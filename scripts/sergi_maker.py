import numpy as np
from os import listdir
from ds import Problem
from qiskit_optimization.translators import from_docplex_mp
from docplex.mp.model_reader import ModelReader
from qiskit_optimization.converters import LinearEqualityToPenalty
import warnings
from qiskit_optimization.algorithms import OptimizationResultStatus

def build_qubo(n_qubs, folder, file, M_strat, towrite_folder):
    m = ModelReader.read(folder + file, ignore_names=True)
    qp = from_docplex_mp(m)
    p = Problem(qp)
    p.qp.name = folder + file

    if M_strat == "qiskit_M":
        converter = LinearEqualityToPenalty()
        nex_folder = "qiskM"
        char = "Q"
    elif M_strat == "our_M":
        M = p.our_M()
        converter = LinearEqualityToPenalty(penalty = M)
        nex_folder = "ourM"
        char = "O"

    qubo = converter.convert(p.qp)
    qubo.name = char + file
    filename = f"../toys/qubos/{towrite_folder}/{nex_folder}/{n_qubs}/{file}"
    qubo.write_to_lp_file(filename)
    #print(f"File {filename} wrote")
    return 

def solve_qubo(n_qubs, folder, file, towrite_folder):
    m = ModelReader.read(folder + file, ignore_names=True)
    qp = from_docplex_mp(m)
    p = Problem(qp)
    p.qp.name = folder + file

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = p.solve_exact()
    if res.status != OptimizationResultStatus.SUCCESS:
        print(f"{filename} results to be infeasible, with status {res.status}")

    fval = np.rint(res.fval)
    x = np.rint(res.x).astype(int)
    filename = f"../toys/qubos/{towrite_folder}/solution/{n_qubs}/{file}"

    f = open(filename, "w")
    f.write(str(x) + "\n")
    f.write(str(fval) + "\n")
    f.close()
    return 


bvars = np.arange(6, 25, 3)
n_samples = 200
M_strategies = ["our_M", "qiskit_M"]
test_set = "../toys/PO_sp500_part3_ra10_mult2"
towrite_folder = "snp_PO"

for i in range(len(bvars)):
    n_qubs = bvars[i]
    print("\n" + str(n_qubs))
    folder = test_set+"/"+str(n_qubs)+"/"
    files = sorted(listdir(folder))
    if len(files) < n_samples:
        raise ValueError(f"Folder {folder} contains only {len(files)} instances, {n_samples} were requested")
    for sample in range(n_samples):
        filename = folder + files[sample]
        print(sample, end = ", ")
        solve_qubo(n_qubs, folder, files[sample], towrite_folder)
        for M_strat in M_strategies:
            build_qubo(n_qubs, folder, files[sample], M_strat, towrite_folder)