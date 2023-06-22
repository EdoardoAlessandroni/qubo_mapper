import numpy as np
import os
from ds import Problem
from qiskit_optimization.translators import from_docplex_mp
from docplex.mp.model_reader import ModelReader
from qiskit_optimization.converters import LinearEqualityToPenalty
import warnings
from qiskit_optimization.algorithms import OptimizationResultStatus

def build_qubos(input_folder, filename, M_strategies, output_folder):
    
    os.mkdir( os.path.join("/home/edo/Desktop/qubo/code/scripts", output_folder+filename[:-3]) )
    os.mkdir( os.path.join("/home/edo/Desktop/qubo/code/scripts", output_folder+filename[:-3]+"/M") )

    m = ModelReader.read(input_folder + filename, ignore_names=True)
    qp = from_docplex_mp(m)
    p = Problem(qp)
    p.qp.name = input_folder + filename

    # write constrained problem
    #qp.name = filename
    os.mkdir( os.path.join("/home/edo/Desktop/qubo/code/scripts", output_folder+filename[:-3]+"/constrained") )
    file = f"../{output_folder}{filename[:-3]}/constrained/{filename}"
    qp.write_to_lp_file(file)

    # write unconstrained QUBOs formulations, with different Ms strategies
    for M_strat in M_strategies:
        os.mkdir( os.path.join("/home/edo/Desktop/qubo/code/scripts", output_folder+filename[:-3]+"/"+M_strat) )
        if M_strat == "qiskit_M":
            converter = LinearEqualityToPenalty()
            char = "Q"
        elif M_strat == "our_M":
            M = p.our_M()
            converter = LinearEqualityToPenalty(penalty = M)
            char = "O"

        qubo = converter.convert(p.qp)
        qubo.name = char + filename
        final_name = filename[:-3]+"_"+M_strat
        file = f"../{output_folder}{filename[:-3]}/{M_strat}/{final_name}"
        qubo.write_to_lp_file(file)
        print(f"File {file} wrote")

        # write M value down
        M = converter.penalty
        f = open(f"../{output_folder}{filename[:-3]}/M/{M_strat}.txt", "w")
        f.write(str(M)+"\n")
        f.close()
    return 


######## MAIN #########

bvars = [100]
n_samples = 20
M_strategies = ["our_M", "qiskit_M"]
#M_strategies = ["our_M"]
#input_test_set = "../../toys/PO_sp500_part3_ra10_mult4"
input_test_set = "../../toys/PO_big"
output_directory = "big_qubos"
#output_directory = "try_dir"

for i in range(len(bvars)):
    n_qubs = bvars[i]
    print("\n" + str(n_qubs))

    output_folder = output_directory+"/"+str(n_qubs)+"/"
    input_folder = input_test_set+"/"+str(n_qubs)+"/"
    filenames = sorted(os.listdir(input_folder))

    for sample, filename in enumerate(filenames[:n_samples]):
        print(sample, end = ", ")
        build_qubos(input_folder, filename, M_strategies, output_folder)

    