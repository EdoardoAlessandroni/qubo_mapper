from os import listdir
import numpy as np
from qiskit_optimization.algorithms import OptimizationResultStatus
from qiskit_optimization.translators import from_docplex_mp
from docplex.mp.model_reader import ModelReader
import pickle
import warnings
from time import time

from ds import Datas, Problem


def run_instance(filename, M_strategies, data, indexes, analyze_gaps):
    '''
    Read LP file to get problem instance, solve it both in a classic and quantum way(s) and compute the gaps
    Return:
        p - the problem instance
        xs - the results got with the M strategies
    '''
    bvars = data.bvars[indexes[0]]
    m = ModelReader.read(filename, ignore_names=True)
    qp = from_docplex_mp(m)
    p = Problem(qp)
    p.qp.name = filename
    data.filenames[indexes[0], indexes[1]] = filename

    # solve classically
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        c_res = p.solve_exact()
    if c_res.status != OptimizationResultStatus.SUCCESS:
        print(f"{filename} results to be infeasible, with status {c_res.status}")
    
    # solve quantumly
    xs = np.ndarray((len(M_strategies), bvars), dtype = int)
    for M_idx in range(len(M_strategies)):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            q_res, M = p.solve_quantum(M_strategies[M_idx])
        data.fval[indexes[0], indexes[1], 0, M_idx] = np.rint(c_res.fval)
        data.fval[indexes[0], indexes[1], 1, M_idx] = np.rint(q_res.fval)
        data.M[indexes[0], indexes[1], M_idx] = M
        xs[M_idx] = np.rint( q_res.x ).astype(int)
        # analyze gaps
        if analyze_gaps:
            H = p.get_obj_hamiltonian()
            Hc = p.get_constraint_hamiltonian()
            H = H + M*Hc # all diagonal
            evs = np.unique(H.round(decimals=14))
            data.gap_norm[indexes[0], indexes[1], M_idx] = (evs[1] - evs[0])/(evs[-1] - evs[0]) # gap of Hamiltonian shifted and squeezed s.t. spectrum is in [0,1]
            data.spec_norm[indexes[0], indexes[1], M_idx] = (evs[-1] - evs[0]) # spectral gap := maximum energy - minimum energy
            data.max_ener[indexes[0], indexes[1], M_idx] = evs[-1]
    return p, xs



def run_test(test_set, bvars, n_samples, M_strategies, analyze_gaps):
    '''
    Run simulation of problems (read from files) for different number of qubits, M-choice strategies, and samples and return data acquired
    '''
    n_M_strategies = len(M_strategies)
    data = Datas(bvars, n_samples, M_strategies)
    # get filenames
    for i in range(len(bvars)):
        n_qubs = bvars[i]
        print(f"\n{n_qubs}")
        folder = f"{test_set}/{n_qubs}/"        
        files = sorted(listdir(folder))
        if len(files) < n_samples:
            raise ValueError(f"Folder {folder} contains only {len(files)} instances, {n_samples} were requested")
        
        # run instances and check feasibility solutions
        tic = time()
        for sample in range(n_samples):
            filename = folder + files[sample]
            print(sample, end = ", ")
            p, xs = run_instance(filename, M_strategies, data, [i, sample], analyze_gaps)
            for bigM_idx in range(n_M_strategies):
                evaluate_feasibility(p, xs[bigM_idx], data, [i, sample, bigM_idx])
        tac = time()
        print(f"It took {tac - tic} sec")
    return data




def run_test_specific_instances(bvars, filenames, n_samples, M_strategies, analyze_gaps):
    '''
    Run simulation of problems (read from files) for different number of qubits, M-choice strategies, and samples and return data acquired. Filenames indicate the filenames of the instances to run, for every bvar
    '''
    n_M_strategies = len(M_strategies)
    data = Datas(bvars, n_samples, M_strategies)

    assert np.shape(filenames) == (len(bvars), n_samples)
    # get filenames
    for i in range(len(bvars)):
        n_qubs = bvars[i]
        print("\n" + str(n_qubs))
        
        # run instances and check feasibility solutions
        tic = time()
        for sample in range(n_samples):
            print(sample, filenames[i][sample])
            p, xs = run_instance(filenames[i][sample], M_strategies, data, [i, sample], analyze_gaps)
            for bigM_idx in range(n_M_strategies):
                evaluate_feasibility(p, xs[bigM_idx], data, [i, sample, bigM_idx])
        tac = time()
        print(f"It took {tac - tic} sec")
    return data


def instances_filenames_biggap(test_set, bvars, n_samples):
    filenames = []
    path = "../toys_adiabevol_simul/" + test_set[15:]
    for i, nvar in enumerate(bvars):
        filenames.append([])
        file_instances = f"{path}/{nvar}/inst_names.txt"

        f = open(file_instances, "r")
        for j in range(n_samples):
            line = f.readline()
            if not line:
                break
            filenames[i].append( line[6:-1] )
        f.close()
    return filenames



def evaluate_feasibility(p, x, data, indexes):
    if p.qp.is_feasible(x): # feasible
        data.is_feasible[indexes[0], indexes[1], indexes[2]] = True
    else: # non feasible
        violated_cons = p.qp.get_feasibility_info(x)[2]
        data.n_violations[indexes[0], indexes[1], indexes[2]] = len(violated_cons)
        max_viol = 0
        for cons in violated_cons:
            violation = np.abs( cons.evaluate(x) - cons.rhs )
            if violation > max_viol:
                max_viol = violation
        data.max_violation[indexes[0], indexes[1], indexes[2]] = max_viol
    return


# run single instance
""" d = Datas([4], 1, 1)
run_instance("../../problems/NN_linear_deg5/4/random10042_4_1.lp", ["our_M"], d, [0,0], True)  """

# ANALYZE DATASET
bvars = np.arange(6, 10, 3)
n_samples = 10
M_strategies = ["our_M", "qiskit_M"]
test_set = "../../problems/NN_linear_deg5"
#test_set = "/home/users/edoardo.alessandroni/codes/problems/PO_norm_part3_mult4"
analyze_gaps = True
data = run_test(test_set, bvars, n_samples, M_strategies, analyze_gaps)



### get big-gap instances filenames
# full_filenames = instances_filenames_biggap(test_set, bvars, n_samples)
# data = run_test_specific_instances(bvars, full_filenames, n_samples, M_strategies, analyze_gaps)


# Save Datas()
#file = open("../../data/NN_linear_deg5_maxener.txt", "wb")
#file = open("/home/users/edoardo.alessandroni/codes/data/test_44.txt", "wb")
#pickle.dump(data, file)
#file.close()