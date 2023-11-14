from os import listdir
import numpy as np
from qiskit_optimization.algorithms import OptimizationResultStatus
from qiskit_optimization.converters import LinearEqualityToPenalty
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
        #data.time[indexes[0], indexes[1], M_idx] = t
        xs[M_idx] = np.rint( q_res.x ).astype(int)
        # analyze gaps
        if analyze_gaps:
            H = p.get_obj_hamiltonian()
            Hc = p.get_constraint_hamiltonian()
            H = H + M*Hc # that's all diagonal!
            evs = np.unique(H.round(decimals=14))   ### THAT WAS CHANGED
            data.gap_norm[indexes[0], indexes[1], M_idx] = (evs[1] - evs[0])/(evs[-1] - evs[0]) # dividing by spectral width gives gap of Hamiltonian shifted and squeezed s.t. spectrum is in [0,1]
    return p, xs


def run_instance_Monly(filename, M_strategies, data, indexes):
    '''
    Read LP file to get problem instance, solve it both in a classic and quantum way(s) and compute the gaps
    Return:
        p - the problem instance
        xs - the results got with the M strategies
    '''
    tic = time()
    bvars = data.bvars[indexes[0]]
    m = ModelReader.read(filename, ignore_names=True)
    tac = time()
    print(f"Model read in {tac - tic} sec")
    qp = from_docplex_mp(m)
    p = Problem(qp)
    p.qp.name = filename
    data.filenames[indexes[0], indexes[1]] = filename
    tic = time()
    print(f"Structure ready in {tic - tac} sec")
    
    # solve quantumly to get M
    for M_idx in range(len(M_strategies)):
        if M_strategies[M_idx] == "heuristic_PO_M":
            M = p.heuristic_PO_M()
        elif M_strategies[M_idx] == "qiskit_M":
            M = np.sum(np.abs(p.obj_quadratic)) + np.sum(np.abs(p.obj_linear)) + 1
        elif M_strategies[M_idx] == "our_M":
            M = p.our_M()

        data.M[indexes[0], indexes[1], M_idx] = M
    return p


# def gaps_various_alphas(n_alphas, H_0, H):
#     alphas = np.linspace(0.9, 1, n_alphas)
#     gaps = np.zeros((n_alphas))
#     for i in range(n_alphas):
#         H_mix = (1-alphas[i])*H_0 + alphas[i]*np.diag(H)
#         evs = np.unique(np.linalg.eigvalsh(H_mix))
#         gaps[i] = (evs[1] - evs[0])/(evs[-1] - evs[0])
#     return gaps
    


def run_test(test_set, bvars, n_samples, M_strategies, analyze_gaps):
    '''
    Run simulation of problems (read from files) for different number of qubits, M-choice strategies, and samples and return data acquired
    '''
    n_M_strategies = len(M_strategies)
    data = Datas(bvars, n_samples, M_strategies)  ### JUST CHANGED TO INTRODUCE M_strategies attribute to Data datastructure
    for i in range(len(bvars)):
        n_qubs = bvars[i]
        print("\n" + str(n_qubs))

        folder = test_set+"/"+str(n_qubs)+"/"
        files = sorted(listdir(folder))
        if len(files) < n_samples:
            raise ValueError(f"Folder {folder} contains only {len(files)} instances, {n_samples} were requested")
        
        tic = time()
        #for sample in range(n_samples):
        #    filename = folder + files[sample]
        #    print(sample, end = ", ")
        #    p, xs = run_instance(filename, M_strategies, data, [i, sample], analyze_gaps)
        #    for bigM_idx in range(n_M_strategies):
        #        evaluate_feasibility(p, xs[bigM_idx], data, [i, sample, bigM_idx])

        for sample in range(n_samples):
            filename = folder + files[sample]
            print(sample, end = ", ")
            p = run_instance_Monly(filename, M_strategies, data, [i, sample])
        tac = time()

        # file_of_filenames = "../easy_toys_adiabevol/"+test_set[11:]+"/easy_inst.txt"
        # filenames = []
        # f = open(file_of_filenames, "r")
        # while True:
        #     line = f.readline()
        #     if not line:
        #         break
        #     filenames.append(line[:-1])
        # f.close()
        # tic = time()
        # for sample, filename in enumerate(filenames):
        #     print(sample, end = ", ")
        #     p, xs = run_instance(filename, M_strategies, data, [i, sample], analyze_gaps)
        #     for bigM_idx in range(n_M_strategies):
        #         evaluate_feasibility(p, xs[bigM_idx], data, [i, sample, bigM_idx])
        # tac = time()

        print(f"It took {tac - tic} sec")
    return data



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


# run single instance
""" d = Datas([4], 1, 1)
run_instance("../../toys/NN_linear_deg5/4/random10042_4_1.lp", ["our_M"], d, [0,0], True)  """


# ANALYZE DATASET
#bvars = np.arange(6, 25, 3)
bvars = [300]
n_samples = 1
#M_strategies = ["our_M", "qiskit_M", "optimal_M", "babbush_M"]
M_strategies = ["heuristic_PO_M", "qiskit_M"]
test_set = "../../toys/PO_big_norm"
analyze_gaps = False
data = run_test(test_set, bvars, n_samples, M_strategies, analyze_gaps)


# Save Datas()
#file = open("../../data/PO_greedy_big_norm_200.txt", "wb")
#pickle.dump(data, file)
#ile.close()