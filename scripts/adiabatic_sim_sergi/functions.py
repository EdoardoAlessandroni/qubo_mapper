import numpy as np
import os
from qibo import hamiltonians
from qibo.symbols import X, Z
from qibo.models.evolution import AdiabaticEvolution
from qibo.models.variational import QAOA
from sympy import symbols
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication_application
from sympy.simplify.simplify import simplify
import matplotlib.pyplot as plt


def get_filenames(folder, typeM, nqubits):
    """Function that collects all the files in a folder into a list for further processing.
    Args:
        folder (str): name of the folder where the instances are saved.
        typeM (str): type of big M computation method. Can be 'solution' to get the solution files.
        nqubits (int): number of qubits to study.
    
    Returns:
        files (list): list with the paths to all the files to study.
    
    """
    files = []
    directory = folder+'/'+typeM+'/'+str(nqubits)
    for filename in os.listdir(directory):
        f = os.path.join(directory, filename)
        files.append(f)
    return files


def get_sym_obj(file, dict_x):
    """Get the symbolic representatio of the instance of a given file.
    Args:
        file (str): path to the file to study.
        dict_x (dict): dictionary needed to parse the expression into symbols.

    Returns:
        obj: sympy expression of the QUBO instance.

    """
    f = open(file, "r")
    obj = []
    for line in f.readlines()[5:]:
        if line == 'Subject To\n':
            break
        obj.append(line[6:-1])
    obj = ' '.join(obj)
    obj = obj.replace('^', '**')
    obj = obj.replace('[', '(')
    obj = obj.replace(']', ')')
    transformations = (standard_transformations + (implicit_multiplication_application,))
    obj = parse_expr(obj, dict_x, transformations=transformations)
    return obj


def substitutions(x):
    """Substitutions of the symbols in the instance into pauli matrices for quantum simulation.
    Args:
        x (list): list of the symbols used in the expression.

    Returns:
        s (list): list of the substitutions needed to go to pauli matrices.

    """
    s = []
    for i in range(len(x)):
        s.append((x[i], (1-Z(i))/2))
    return s


def get_maxmin_eigenvalues(expression):
    """Compute the minimum and maximum eigenvalues for a given expression for normalization.
    Args:
        expression: sympy expression to compute the minimum and maximum eigenvalues.

    Returns:
        min, max (float): minimum and maximum eigenvalues of the expression.

    """
    m = hamiltonians.SymbolicHamiltonian(expression)
    m_eigs = m.eigenvalues()
    mx = max(m_eigs)
    mn = min(m_eigs)
    return mx, mn


def get_sym_ham(sym_obj, x):
    """Get symbolic hamiltonian from sympy expression, normalized.
    Args:
        sym_obj: symbolic expression of the QUBO instance.
        x (list): list of symbols that appear in the expression

    Returns:
        ham (SymbolicHamiltonian): hamiltonian qibo object that encodes the qubo instance.

    """
    subs = substitutions(x)
    ham = simplify(sym_obj.subs(subs))
    max_eig, min_eig = get_maxmin_eigenvalues(ham)
    ham = hamiltonians.SymbolicHamiltonian(simplify((ham-min_eig)/(max_eig-min_eig)))
    return ham


def get_h0(nqubits):
    """Generate the initial hamiltonian for the adiabatic evolution/qaoa, normalized.
    Args:
        nqubits (int): number of qubits of the instance.

    Returns:
        ham0 (SymbolicHamiltonian): X hamiltonian used for initial hamiltonian or mixer hamiltonian.

    """
    ham0 = sum((0.5 * (1 - X(i))) for i in range(nqubits))
    ham0 = hamiltonians.SymbolicHamiltonian(simplify(ham0/nqubits))
    return ham0


def get_solution(file, nqubits):
    """Get the solution of the instance from the solution files.
    Args:
        file (str): path to the solution file to analyze.
        nqubits (int): number of qubits of the instance to go from binary to int.

    Returns:
        sol (int): number of the solution as if it was in binary.
        sol_bin (str): binary configuration of the solution.

    """
    f = open(file, "r")
    sol_bin = f.readline()
    sol_bin = sol_bin[1:-1:2]
    sol = 0
    for i in range(len(sol_bin)):
        sol += int(sol_bin[i])*2**(nqubits-1-i)

    return sol, sol_bin


def compare_adiabatic(files_our, files_qiskit, files_solution, T, nqubits, dt=.1, solver='exp'):
    """Funtion that compares ourM and qiskitM formulation using adiabatic evolution for a fixed final time.
    Args:
        files_our (list): list of files with ourM to be studied.
        files_qiskit (list): list of files with qiskitM to be studied.
        files_solution (list): list of files with the solution of the instances.
        T (int): final time for the adiabatic evolution.
        nqubits (int): number of qubits of the instance.
        dt (float): step for the integration of the adiabatic path.
        solver (str): solver to use for the simulation.

    Returns:
        probs_ours (np.array): array of probabilities of measuring the target state with ourM.
        probs_qiskit (np.array): array of probabilities of measuring the target state with qiskitM.

    """
    x = symbols(" ".join((f"x{i}" for i in range(0, nqubits))))
    dict_x = {str(xx):xx for xx in x}

    probs_ours = []
    probs_qiskit = []

    s = lambda t: t

    for i in range(len(files_our)):
        sym_obj = get_sym_obj(files_our[i], dict_x)
        sym_ham = get_sym_ham(sym_obj, x)
        ham0 = get_h0(nqubits)
        sol, sol_bin = get_solution(files_solution[i], nqubits)
        evolve = AdiabaticEvolution(ham0, sym_ham, s, dt=dt, solver=solver)
        final_state = evolve(final_time=T)

        probs_ours.append((abs(final_state[sol])**2))

        sym_obj = get_sym_obj(files_qiskit[i], dict_x)
        sym_ham = get_sym_ham(sym_obj, x)
        ham0 = get_h0(nqubits)
        evolve = AdiabaticEvolution(ham0, sym_ham, s, dt=dt, solver=solver)
        final_state = evolve(final_time=T)
        
        probs_qiskit.append((abs(final_state[sol])**2))
    return np.array(probs_ours), np.array(probs_qiskit)


def compare_qaoa(files_our, files_qiskit, files_solution, nqubits, nparams):
    """Funtion that compares ourM and qiskitM formulation using adiabatic evolution for a fixed final time.
    Args:
        files_our (list): list of files with ourM to be studied.
        files_qiskit (list): list of files with qiskitM to be studied.
        files_solution (list): list of files with the solution of the instances.
        nqubits (int): number of qubits of the instance.
        nparams (str): number of parameters for the optimization. Will always be 2*nparams as it has to be even.

    Returns:
        probs_ours (np.array): array of probabilities of measuring the target state with ourM.
        probs_qiskit (np.array): array of probabilities of measuring the target state with qiskitM.

    """
    x = symbols(" ".join((f"x{i}" for i in range(0, nqubits))))
    dict_x = {str(xx):xx for xx in x}

    ham0 = get_h0(nqubits)

    probs_ours = []
    probs_qiskit = []
    time_ours = []
    time_qiskit = []

    s = lambda t: t

    initial_parameters = 0.01 * np.random.random(2*nparams)

    for i in range(len(files_our)):
        sym_obj = get_sym_obj(files_our[i], dict_x)
        sym_ham = get_sym_ham(sym_obj, x)
        sol, sol_bin = get_solution(files_solution[i], nqubits)

        qaoa = QAOA(sym_ham, mixer=ham0, solver='exp')
        best_energy, final_parameters, extra = qaoa.minimize(initial_parameters,  
                                                             method="L-BFGS-B")
        
        qaoa.set_parameters(final_parameters)
        final_state = qaoa.execute()

        probs_ours.append((abs(final_state[sol])**2))
        
        sym_obj = get_sym_obj(files_qiskit[i], dict_x)
        sym_ham = get_sym_ham(sym_obj, x)

        qaoa = QAOA(sym_ham, mixer=ham0, solver='exp')
        best_energy, final_parameters, extra = qaoa.minimize(initial_parameters, 
                                                             method="L-BFGS-B")
        
        qaoa.set_parameters(final_parameters)
        final_state = qaoa.execute()

        probs_qiskit.append((abs(final_state[sol])**2))

    return np.array(probs_ours), np.array(probs_qiskit)


def save_data(data, name_data):
    """General function to save some data.
    Args:
        data: data to be saved.
        name_data (str): name for the file where the data is to be saved.

    """
    np.save(f'data/{name_data}.npy', data)
