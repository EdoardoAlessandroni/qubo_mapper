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
        M_folder = "qiskit_M"
        char = "Q"
    elif M_strat == "our_M":
        M = p.our_M()
        converter = LinearEqualityToPenalty(penalty = M)
        M_folder = "our_M"
        char = "O"
    elif M_strat == "babbush_M":
        M = p.babbush_M()
        converter = LinearEqualityToPenalty(penalty = M)
        M_folder = "babbushM"
        char = "B"

    qubo = converter.convert(p.qp)
    qubo.name = char + file
    filename = f"../{towrite_folder}/{M_folder}/{n_qubs}/{file}"
    qubo.write_to_lp_file(filename)
    print(f"File {filename} wrote")
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

    filename = f"../{towrite_folder}/solution/{n_qubs}/{file}"

    f = open(filename, "w")
    f.write(str(x) + "\n")
    f.write(str(fval) + "\n")
    f.close()

    # write constrained problem as well
    qp.name = file
    filename = f"../{towrite_folder}/constrained/{n_qubs}/{file}"
    qp.write_to_lp_file(filename)
    return 


def post_process_filenames(filenames):
    i = 0
    while i < len(filenames[0]):
        if filenames[0][i:i+4] != "rand":
            i += 1
        else:
            break
    new_filenames = []
    for file in filenames:
        new_filenames.append(file[i:])
    return new_filenames






######## MAIN #########

bvars = np.arange(21, 22)
n_samples = 10
M_strategies = ["our_M", "qiskit_M"]
test_set = "../../problems/PO_sp500_part3_ra10_mult2"
towrite_folder = "easy_toys_adiabevol/PO_sp500_part3_ra10_mult2"

for i in range(len(bvars)):
    n_qubs = bvars[i]
    print("\n" + str(n_qubs))

    folder = test_set+"/"+str(n_qubs)+"/"
    file_of_filenames = "../easy_toys_adiabevol/"+test_set[11:]+"/easy_inst21.txt"
    filenames = []
    f = open(file_of_filenames, "r")
    while True:
        line = f.readline()
        if not line:
            break
        filenames.append(line[:-1])
    f.close()
    filenames = sorted(post_process_filenames(filenames))

    for sample, filename in enumerate(filenames[:n_samples]):
        print(sample, end = ", ")
        solve_qubo(n_qubs, folder, filename, towrite_folder)
        for M_strat in M_strategies:
            build_qubo(n_qubs, folder, filename, M_strat, towrite_folder)