import numpy as np
from os import listdir
from ds import Problem
from qiskit_optimization.translators import from_docplex_mp
from docplex.mp.model_reader import ModelReader
from qiskit_optimization.converters import LinearEqualityToPenalty
import warnings
from qiskit_optimization.algorithms import OptimizationResultStatus
import pickle

def build_qubo(n_qubs, folder, file, M_strat, towrite_folder):
    ''' Builds the QUBO with a given Big-M strategy and writes it on file '''
    m = ModelReader.read(folder + file, ignore_names=True)
    qp = from_docplex_mp(m)
    p = Problem(qp)
    p.qp.name = folder + file

    if M_strat == "qiskitM":
        converter = LinearEqualityToPenalty()
        char = "Q"
    elif M_strat == "ourM":
        M = p.our_M()
        converter = LinearEqualityToPenalty(penalty = M)
        char = "O"
    elif M_strat == "babbushM":
        M = p.babbush_M()
        converter = LinearEqualityToPenalty(penalty = M)
        char = "B"

    # write down QUBO reformulation
    qubo = converter.convert(p.qp)
    qubo.name = char + file
    filename = f"../{towrite_folder}/{n_qubs}/{M_strat}/QUBO/{file}"
    qubo.write_to_lp_file(filename)

    # write down spectral gap (and maximal energy)
    gap, max_ener = spectral_gap(n_qubs, file, M_strat)
    filename = f"../{towrite_folder}/{n_qubs}/{M_strat}/maxener/{file}"
    f = open(filename, "w")
    f.write(str(max_ener)+"\n")
    f.close()
    filename = f"../{towrite_folder}/{n_qubs}/{M_strat}/gap/{file}"
    f = open(filename, "w")
    f.write(str(gap)+"\n")
    f.close()

    print(f"File {filename} and its QUBO reform. wrote")
    return 

def solve_qubo(n_qubs, folder, file, towrite_folder):
    ''' Writes the 'constrained' problem on file, along with its solution (x, f(x)), in another file '''
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

    filename = f"../{towrite_folder}/{n_qubs}/solution/{file}"
    f = open(filename, "w")
    f.write(f"x = {x}\nf(x) = {fval}\n")
    f.close()

    # write constrained problem as well
    qp.name = file
    filename = f"../{towrite_folder}/{n_qubs}/constrained/{file}"
    qp.write_to_lp_file(filename)
    print(f"Constrained and solution wrote")
    return 


def spectral_gap(n_qubs, file, M_strat):
    datafile = open("../../data/" + problem_set + "_maxener.txt", "rb")       ############ "_maxener" is added so the Data instance open and analyzed is the one that not only has the spectral gap, but the maximum energy (and spectral norm) too
    data = pickle.load(datafile)
    datafile.close()
    
    # # set var_idx
    # if problem_set == "NN_linear_deg5":
    #     var_idx = n_qubs - 4
    # elif problem_set == "SPP_p25":
    #     var_idx = n_qubs - 6
    # elif problem_set == "PO_norm_part3_mult4":
    #     var_idx = np.rint( (n_qubs - 6)/3 ).astype(int)
    var_idx = int(n_qubs/6 - 1)

    # set M_strat
    if M_strat == "ourM":
        M_idx = 0
    elif M_strat == "qiskitM":
        M_idx = 1

    # set inst_idx
    tail_len = 13
    for inst_idx, name in enumerate(data.filenames[var_idx]):
        if name[-tail_len:] == file[-tail_len:]:
            #print("Match:", name[-tail_len:], file[-tail_len:])
            break
        if inst_idx == len(data.filenames[var_idx]) - 1:
            print(f"No filename {file[-tail_len:]} found to match")
            raise ValueError()

    gap = data.gap_norm[var_idx, inst_idx, M_idx]
    max_ener = data.max_ener[var_idx, inst_idx, M_idx]
    return gap, max_ener


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

bvars = np.arange(6, 25, 6)
#bvars = np.arange(12, 13, 6)
n_samples = 25
M_strategies = ["ourM", "qiskitM"]

problem_set = "NN_linear_deg5"
towrite_folder = "toys_adiabevol_simul/" + problem_set



for i in range(len(bvars)):
    n_qubs = bvars[i]
    print("\n\n\n" + str(n_qubs))

    folder_instances = f"../../problems/{problem_set}/{n_qubs}/"
    file_of_filenames = f"../toys_adiabevol_simul/{problem_set}/{n_qubs}/inst_names.txt"
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
        print(f"{sample}:")
        solve_qubo(n_qubs, folder_instances, filename, towrite_folder)
        for M_strat in M_strategies:
            build_qubo(n_qubs, folder_instances, filename, M_strat, towrite_folder)






###### Note: 
# to run on Data when only specific number of bvars where run and saved, you should change the var_idx index in file spectral_gap.
# In particular, for data like PO_norm_part3_mult4_maxener.txt, only [6, 12, 18, 24] vars were analyzed, instead in PO_norm_part3_mult4.txt all vars in range(6, 25, 3) were analyzed,
# so adapt this var_idx to the bvars you're analyzing here, together with the bvars run and saved in Data instances