from os import listdir
import numpy as np
from qiskit_optimization.algorithms import OptimizationResultStatus, CplexOptimizer
from qiskit_optimization.translators import from_docplex_mp
from qiskit_optimization.converters import LinearEqualityToPenalty
from docplex.mp.model_reader import ModelReader
from pickle import dump as Pdump
import warnings
from time import time
from copy import deepcopy

from ds_ame import Datas, Problem, Converter




def run_instance(filename, data, indexes, analyze_gaps):
    '''
    Read LP file to get problem instance, solve it iteratively in the AME approach
    Return:
        p - the problem instance
    '''
    mod = ModelReader.read(filename, ignore_names=True)
    qp = from_docplex_mp(mod)
    p = Problem(qp)
    p.qp.name = filename
    m = p.n_const

    # solve classically
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        c_res = p.solve_exact()
    if c_res.status != OptimizationResultStatus.SUCCESS:
        #print(f"{filename} results to be infeasible!")
        print(f"{filename} results to be infeasible, with status {c_res.status}")
    data.fval_classic[indexes[0], indexes[1]] = np.rint(c_res.fval).astype(int)


    # solve iteratively with AME
    converter = Converter(p)
    our_M = p.our_M()
    initialM_divider = 10
    M = np.rint(our_M / initialM_divider).astype(int)
    if M == 0:
        M = 1
    M_ls = M*np.ones(m)   
    is_feas = False
    fvals, gaps = [], []
    viol_nums, Ms = [],[] # list [len=algo_steps] of arrays [len=num_constr]
    step = 0
    while not is_feas:
        # solve Qubo
        qubo = converter.convert(M_ls)
        res = CplexOptimizer(disp = False).solve(qubo)

        # check solution
        f, x = np.rint(res.fval).astype(int), np.rint( res.x ).astype(int)
        is_feas = p.qp.is_feasible(x)

        # get violation numbers
        k = np.zeros(m)
        cons = p.qp.linear_constraints
        for i, c in enumerate(cons):
            k[i] = (c.evaluate(x) - c.rhs )**2
            k[i] = np.rint(k[i]).astype(int)

        # save info
        viol_nums.append(deepcopy(k))
        Ms.append(deepcopy(M_ls))
        fvals.append(deepcopy(f))
        if analyze_gaps:
            p_qubo = Problem(qubo)
            H = p_qubo.get_obj_hamiltonian()
            evs = np.unique(H) # already sorted
            gap = (evs[1] - evs[0]) / (evs[-1] - evs[0])
            gaps.append(gap) # dividing by spectral width gives gap of Hamiltonian shifted and squeezed s.t. spectrum is in [0,1]

        # map to next step
        for j in range(m):
            M_ls[j] = min((k[j] + 1)*M_ls[j], our_M*1.1)
        step += 1

    our_qubo = converter.convert( our_M*np.ones(m) )
    our_qubo = Problem(our_qubo)
    H = our_qubo.get_obj_hamiltonian()
    evs = np.unique(H) # already sorted
    our_gap = (evs[1] - evs[0]) / (evs[-1] - evs[0])

    # fill data with relevant info
    data.num_const[indexes[0],indexes[1]] = p.n_const
    data.max_iter[indexes[0],indexes[1]] = step
    data.gaps[f"{indexes[0]}_{indexes[1]}"] = gaps
    data.our_gap[indexes[0],indexes[1]] = our_gap
    data.fvals[f"{indexes[0]}_{indexes[1]}"] = np.array((fvals))
    data.Ms[f"{indexes[0]}_{indexes[1]}"] = Ms
    data.violation_nums[f"{indexes[0]}_{indexes[1]}"] = viol_nums
    data.is_optimum[indexes[0],indexes[1]] = ((data.fval_classic[indexes[0], indexes[1]] - fvals[-1]) == 0 )
    return p
    


def run_test(test_set, bvars, n_samples, analyze_gaps):
    '''
    Run simulation of problems (read from files) for different number of qubits, M-choice strategies, and samples and return data acquired
    '''
    data = Datas(bvars, n_samples)
    for i in range(len(bvars)):
        n_qubs = bvars[i]
        print("\n" + str(n_qubs))
        folder = test_set+"/"+str(n_qubs)+"/"
        files = sorted(listdir(folder))
        if len(files) < n_samples:
            raise ValueError(f"Folder {folder} contains only {len(files)} instances, {n_samples} were requested")
        
        tic = time()
        for sample in range(n_samples):
            filename = folder + files[sample]
            print(sample, end = ", ")
            p = run_instance(filename, data, [i, sample], analyze_gaps)
        tac = time()
        m_max = p.n_const
        print(f"It took {np.round(tac - tic, 4)} sec")

    # postprocess dictionaries into big arrays
    create_arr_from_dicts(data, m_max)
    return data

def create_arr_from_dicts(data, m_max):
    """ Now that max_number_step_ame is known, create arrays to store relevent data from unstructured dictionaries """
    max_iter_tot = data.max_iter.max().astype(int)
    n_bvars, n_samples = np.shape(data.max_iter)
    
    Ms, fvals = deepcopy(data.Ms), deepcopy(data.fvals)
    gaps, ks = deepcopy(data.gaps), deepcopy(data.violation_nums)
    data.fvals, data.gaps = np.zeros((n_bvars, n_samples, max_iter_tot)), np.zeros((n_bvars, n_samples, max_iter_tot))
    data.Ms, data.violation_nums = np.zeros((n_bvars, n_samples, max_iter_tot, m_max)), np.zeros((n_bvars, n_samples, max_iter_tot, m_max))
    
    for i_vars in np.arange(n_bvars):
        for i_sample in np.arange(n_samples):
            key = f"{i_vars}_{i_sample}"
            iter_idx = data.max_iter[i_vars, i_sample]
            m = data.num_const[i_vars, i_sample]
            data.fvals[i_vars, i_sample, :iter_idx] = fvals[key]
            data.gaps[i_vars, i_sample, :iter_idx] = gaps[key]
            for t in range(iter_idx):
                data.Ms[i_vars, i_sample, t, :m] = Ms[key][t]
                data.violation_nums[i_vars, i_sample, t, :m] = ks[key][t]



# run single instance
""" d = Datas([10], 1, 2)
run_instance(filename, data, indexes, analyze_gaps) """


# ANALYZE DATABASE
bvars = np.arange(20, 26)
n_samples = 200
test_set = "../../toys/SPP_p15"
analyze_gaps = True
data = run_test(test_set, bvars, n_samples, analyze_gaps)


# Save Datas()
file = open("../../data/ame_manym/SPP_p15_from20.txt", "wb")
Pdump(data, file)
file.close()