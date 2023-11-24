from os import listdir
import numpy as np
from qiskit_optimization.algorithms import OptimizationResultStatus
from qiskit_optimization.converters import LinearEqualityToPenalty
from qiskit_optimization.translators import from_docplex_mp
from docplex.mp.model_reader import ModelReader
import pickle
from time import time

import multiprocessing as mp
from itertools import product

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



def run_test(test_set, bvars, n_samples, M_strategies):
    '''
    Run simulation of problems (read from files) for different number of qubits, M-choice strategies, and samples and return data acquired
    '''
    data = Datas(bvars, n_samples, M_strategies)
    for i in range(len(bvars)):
        n_qubs = bvars[i]
        print("\n" + str(n_qubs))

        folder = test_set+"/"+str(n_qubs)+"/"
        files = sorted(listdir(folder))
        if len(files) < n_samples:
            raise ValueError(f"Folder {folder} contains only {len(files)} instances, {n_samples} were requested")
        filenames = [folder + fl for fl in files[:n_samples]]

        args_iterable = product(zip(filenames, range(n_samples)), [M_strategies], [data], [i])

        tic = time()
        with mp.Pool() as pool:
            results = pool.starmap(analyze_instance, args_iterable)
        tac = time()
        print(f"It took {tac - tic} sec")
    return data

def analyze_instance(file_n_sample, M_strategies, data, i):
    filename, sample = file_n_sample
    print(sample, end = ", ")
    return run_instance_Monly(filename, M_strategies, data, [i, sample])

# ANALYZE DATASET
#bvars = np.arange(6, 25, 3)
bvars = [30]
n_samples = 10
M_strategies = ["heuristic_PO_M", "qiskit_M"]
test_set = "/home/users/edoardo.alessandroni/toys/PO_big_norm"
data = run_test(test_set, bvars, n_samples, M_strategies)


# Save Datas()
#file = open("/home/users/edoardo.alessandroni/codes/data/PO_greedy_big_norm_300_second.txt", "wb")
#pickle.dump(data, file)
#file.close()