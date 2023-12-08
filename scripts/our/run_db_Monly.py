from os import listdir
import numpy as np
from qiskit_optimization.algorithms import OptimizationResultStatus
from qiskit_optimization.converters import LinearEqualityToPenalty
from qiskit_optimization.translators import from_docplex_mp
from docplex.mp.model_reader import ModelReader
import pickle
from time import time

from ds import Datas, Problem



def run_instance_Monly(filename, M_strategies, data, indexes):
    '''
    Read LP file to get problem instance and compute the corresponding M only
    Return:
        p - the problem instance
    '''
    m = ModelReader.read(filename, ignore_names=True)
    qp = from_docplex_mp(m)
    p = Problem(qp)
    p.qp.name = filename
    data.filenames[indexes[0], indexes[1]] = filename
    
    # compute M
    for M_idx in range(len(M_strategies)):
        if M_strategies[M_idx] == "heuristic_PO_M":
            M = p.heuristic_PO_M()
        elif M_strategies[M_idx] == "qiskit_M":
            M = np.sum(np.abs(p.obj_quadratic)) + np.sum(np.abs(p.obj_linear)) + 1
        elif M_strategies[M_idx] == "our_M":
            M = p.our_M()

        data.M[indexes[0], indexes[1], M_idx] = M
    return p



def run_test(test_set, bvars, n_samples, M_strategies):
    '''
    Run simulation of problems (read from files) for different number of qubits, M-choice strategies, and samples and return data acquired
    '''
    data = Datas(bvars, n_samples, M_strategies)
    # get filenames
    for i in range(len(bvars)):
        n_qubs = bvars[i]
        print("\n" + str(n_qubs))
        folder = test_set+"/"+str(n_qubs)+"/"
        files = sorted(listdir(folder))
        if len(files) < n_samples:
            raise ValueError(f"Folder {folder} contains only {len(files)} instances, {n_samples} were requested")
        
        # run instances
        tic = time()
        for sample in range(n_samples):
            filename = folder + files[sample]
            print(sample, end = ", ")
            p = run_instance_Monly(filename, M_strategies, data, [i, sample])
        tac = time()
        print(f"It took {tac - tic} sec")
    return data


# run single instance
""" d = Datas([4], 1, 1)
run_instance("../../toys/NN_linear_deg5/4/random10042_4_1.lp", ["our_M"], d, [0,0], True)  """


# ANALYZE DATASET
#bvars = np.arange(6, 25, 3)
bvars = [10]
n_samples = 5
M_strategies = ["heuristic_PO_M"]
test_set = "../../toys/PO_big_norm"
#test_set = "/home/users/edoardo.alessandroni/codes/toys/PO_big_norm"
data = run_test(test_set, bvars, n_samples, M_strategies)


# Save Datas()
#file = open("../../data/test_sdp_big.txt", "wb")
#file = open("/home/users/edoardo.alessandroni/codes/data/PO_greedy_big_norm_75.txt", "wb")
#pickle.dump(data, file)
#file.close()