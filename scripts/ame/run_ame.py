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

from ds_ame import Datas, Problem




def run_instance(filename, data, indexes, analyze_gaps):
    '''
    Read LP file to get problem instance, solve it iteratively in the AME approach
    Return:
        p - the problem instance
    '''
    m = ModelReader.read(filename, ignore_names=True)
    qp = from_docplex_mp(m)
    p = Problem(qp)
    p.qp.name = filename

    if analyze_gaps:
        H = p.get_obj_hamiltonian()
        Hc = p.get_constraint_hamiltonian()

    # solve classically
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        c_res = p.solve_exact()
    if c_res.status != OptimizationResultStatus.SUCCESS:
        #print(f"{filename} results to be infeasible!")
        print(f"{filename} results to be infeasible, with status {c_res.status}")
    data.fval_classic[indexes[0], indexes[1]] = np.rint(c_res.fval).astype(int)

    # solve iteratively with AME
    our_M = p.our_M()
    initialM_divider = 10
    M = our_M / initialM_divider
    #M = np.rint(our_M / initialM_divider).astype(int)
    #if M == 0:
    #    M = 1
    is_feas = False
    fvals, viol_nums, Ms, gaps = [],[],[],[]
    step = 0
    while not is_feas:
        # solve Qubo
        converter = LinearEqualityToPenalty(penalty = M)
        qubo = converter.convert(p.qp)
        res = CplexOptimizer(disp = False).solve(qubo)
        # check solution
        f, x = np.rint(res.fval).astype(int), np.rint( res.x ).astype(int)
        is_feas = p.qp.is_feasible(x)
        # get violation number
        k = 0
        violated_cons = p.qp.get_feasibility_info(x)[2]
        for cons in violated_cons:
            k += ( cons.evaluate(x) - cons.rhs )**2
        k = np.rint(k).astype(int)
        # save info
        viol_nums.append(k)
        Ms.append(M)
        fvals.append(f)
        if analyze_gaps:
            H_tot = H + M*Hc
            evs = np.unique(H_tot.round(decimals=14))   ### THAT WAS CHANGED
            #evs = np.unique(H) # already sorted
            gaps.append((evs[1] - evs[0]) / (evs[-1] - evs[0])) # dividing by spectral width gives gap of Hamiltonian shifted and squeezed s.t. spectrum is in [0,1]
        # map to next step
        M = min((k + 1)*M, our_M)
        step += 1
    
    # compute our_gap
    if analyze_gaps:
        H_tot = H + our_M*Hc
        evs = np.unique(H_tot.round(decimals=14))
        our_gap = (evs[1] - evs[0]) / (evs[-1] - evs[0])
        
    # fill data with relevant info
    data.max_iter[indexes[0],indexes[1]] = step
    data.gaps[f"{indexes[0]}_{indexes[1]}"] = gaps
    data.our_M[indexes[0],indexes[1]] = our_M
    data.our_gap[indexes[0],indexes[1]] = our_gap
    data.Ms[f"{indexes[0]}_{indexes[1]}"] = np.array((Ms)) # or = Ms
    data.fvals[f"{indexes[0]}_{indexes[1]}"] = np.array((fvals))
    data.violation_nums[f"{indexes[0]}_{indexes[1]}"] = np.array((viol_nums))
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
        print(f"It took {tac - tic} sec")

    # postprocess dictionaries into big arrays
    create_arr_from_dicts(data)
    return data

def create_arr_from_dicts(data):
    """ Now that max_number_step_ame is known, create arrays to store relevent data from unstructured dictionaries """
    max_iter_tot = data.max_iter.max().astype(int)
    n_bvars, n_samples = np.shape(data.max_iter)
    
    Ms, fvals = deepcopy(data.Ms), deepcopy(data.fvals)
    gaps, ks = deepcopy(data.gaps), deepcopy(data.violation_nums)
    data.Ms, data.fvals = np.zeros((n_bvars, n_samples, max_iter_tot)), np.zeros((n_bvars, n_samples, max_iter_tot))
    data.gaps, data.violation_nums = np.zeros((n_bvars, n_samples, max_iter_tot)), np.zeros((n_bvars, n_samples, max_iter_tot))
    
    for i_vars in np.arange(n_bvars):
        for i_sample in np.arange(n_samples):
            key = f"{i_vars}_{i_sample}"
            iter_idx = data.max_iter[i_vars, i_sample]
            data.Ms[i_vars, i_sample, :iter_idx] = Ms[key]
            data.fvals[i_vars, i_sample, :iter_idx] = fvals[key]
            data.gaps[i_vars, i_sample, :iter_idx] = gaps[key]
            data.violation_nums[i_vars, i_sample, :iter_idx] = ks[key]



# run single instance
""" d = Datas([10], 1, 2)
run_instance("../toys/NN_25/10/random1042_10_2.lp", ["qiskit_M", "our_M"], d, [0,0], True, True, True)  """


# ANALYZE DATABASE
bvars = np.arange(25, 26)
n_samples = 200
test_set = "../../toys/SPP_p15"
analyze_gaps = True
data = run_test(test_set, bvars, n_samples, analyze_gaps)

# Save Datas()
file = open("../../data/ame/SPP_p15_newM0_25.txt", "wb")
Pdump(data, file)
file.close()