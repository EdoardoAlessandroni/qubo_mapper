import numpy as np
from qiskit_optimization import QuadraticProgram
from qiskit_optimization.algorithms import OptimizationResultStatus
from docplex.mp.model_reader import ModelReader
from qiskit_optimization.translators import from_docplex_mp
import warnings
from os import listdir
from os import remove as osremove

from ds import Problem


# generate random QuadraticProgram
def build_qp(n : int, h, l = None, A = None, b = None) -> QuadraticProgram:
    '''
    Returns a quadratic problem with n decision variable, h and l respectively quadratic and linear 
    part of objective function and A and b matrix and vector corresponding to constraints
    '''

    qp = QuadraticProgram(name='random_problem.lp')
    for i in range(n):
        qp.binary_var()

    # objective function
    qp.minimize(linear = l, quadratic = h)

    # constraints
    if not (A is None):
        n_constr = np.shape(A)[0]
        for i in range(n_constr):
            qp.linear_constraint(linear = A[i], sense = "==", rhs = b[i])
    return qp

""" def rand_sym_matrix_int_old(size): # it was used in BN normalization (in generating) and gap analysis (in running)
    m = np.random.randint(low = -10, high = 11, size = (size, size))
    return m + m.transpose() """

def rand_sym_matrix_density(size, density):
    if density is None:
        raise ValueError("Density has not been given")
    m = np.triu( np.random.randint(low = -10, high = 11, size = (size, size)) )
    indexes = [(i,j) for i in range(size) for j in range(i, size)]
    np.random.shuffle(indexes)
    n_indexes = int( len(indexes)*(1-density) )
    for ij in range(n_indexes):
        m[indexes[ij]] = 0
    return m

def rand_vector_density(size, density):
    if density is None:
        raise ValueError("Density has not been given")
    L = np.random.randint(low = -10, high = 11, size = (size))
    vars_zero = np.random.choice(np.arange(size), size = int((1-density)*size), replace = False)
    L[vars_zero] = 0
    return L

def rand_sym_matrix_degree(size, max_degree):
    if max_degree is None:
        raise ValueError("degree has not been given")
    if max_degree >= size:
        max_degree = size-1

    # algorithm t
    # o generate a random quasi-regular graph of max_degree
    zeroset = [(i,j) for i in range(size) for j in range(i+1,size)]
    rem_links = (max_degree*np.ones(size)).astype("int")
    for i in range(size-1):
        p = [rem_links[j] for j in range(i+1, size)]
        if rem_links[i] == 0 or np.sum(p) == 0:
            continue
        p = p/np.sum(p)
        neighbors = np.random.choice([j for j in range(i+1, size)], size = min(rem_links[i], np.count_nonzero(p)), replace = False, p = p)
        rem_links[neighbors] += -1
        for j in neighbors:
            zeroset.remove((i,j))

    coeffs = np.append(np.arange(-10,0), np.arange(1,11))
    m = np.triu(  np.random.choice(coeffs,  size = (size, size)), k=1)
    for ij in zeroset:
        m[ij] = 0
    return m


""" def rand_sym_matrix_SN(size):
    k = 4/3/(size**2 + size)
    return np.random.uniform(low = -k, high = k, size = (size, size)) 

def rand_vector_SN(size):
    k = 4/3/(size**2 + size)
    return np.random.uniform(low = -k, high = k, size = (size)) """

def rand_constr_matrix_density(size, n_constr, density):
    if density is None:
        density = .5
    A = np.zeros((n_constr, size), dtype=int)
    max = 2
    coef_list = np.delete(np.arange(-max,max+1), obj = max)
    for i in range(n_constr):
        vars_nz = np.random.choice(np.arange(size), size = int(density*size), replace = False)
        coefs = np.random.choice(coef_list, size = int(density*size))
        A[i, vars_nz] = coefs
    return A

def rand_constr_matrix_degree(size, n_constr, max_degree):
    if max_degree is None:
        raise ValueError("degree has not been given")
    if max_degree >= size:
        max_degree = size-1

    density = max_degree/size # density simply given by degree (always <1)
    A = np.zeros((n_constr, size), dtype=int)
    max = 2
    coef_list = np.delete(np.arange(-max,max+1), obj = max)
    for i in range(n_constr):
        vars_nz = np.random.choice(np.arange(size), size = int(density*size), replace = False)
        coefs = np.random.choice(coef_list, size = int(density*size))
        A[i, vars_nz] = coefs
    return A


def rand_constr_matrix_SPP(size, n_constr, max_degree, p_baseline = 0.15):
    """ if max_degree is None:
        raise ValueError("degree has not been given") """
    # hardcoded parameter
    prob = p_baseline
    """ if 19 < size: # hardcoded number which works only for delts = 5 and p_base = 0.35!!
        # inverted formula which statistically enforce maximum desired degree on A^tA which appears on penalized objective
        prob = np.sqrt( 1 - np.power( 1-max_degree/(size-1) , 1/n_constr) ) 
    """
    A = np.zeros((n_constr, size), dtype=int)
    for i in range(size):
        column = np.random.choice([0,1], size = n_constr, replace = True, p = [1-prob, prob])
        # ensure that every column is not trivially 0, which corresponds to the empty subset as a useless variable
        while np.sum(column) == 0: 
            column = np.random.choice([0,1], size = n_constr, replace = True, p = [1-prob, prob])
        A[:, i] = column
    
    # check that we didn't add trivially infeasible constraint of the form 0=1, ie every row must have at least a one
    for i in range(n_constr):
        if np.sum(A[i]) == 0:
            index = np.random.randint(0, size)
            A[i, index] = 1
    return A


def rand_constr_vector(n_constr):
    return np.random.randint(low = -2, high = 2, size = (n_constr))

def random_qp(n_vars, constr : bool, n_constr, density, max_degree, sampling_normalization = False, SPP = False):
    '''
    Define a random QuadraticProgram, with n_vars number of binary variables and n_constr
    number of constraints if given, a random number in [0, n_vars/3] otherwise
    '''
    """ if sampling_normalization:
        h = rand_sym_matrix_SN(n_vars)
        l = rand_vector_SN(n_vars)
    else: """
    if SPP:
        h = np.zeros((n_vars, n_vars))
        l = np.random.randint(low = 1, high = 99, size = (n_vars))
    elif density is not None:
        h = rand_sym_matrix_density(n_vars, density)
        l = rand_vector_density(n_vars, density)
    elif max_degree is not None:
        h = rand_sym_matrix_degree(n_vars, max_degree)
        l = rand_vector_density(n_vars, density = 1)
    else:
        raise ValueError("Both 'density' and 'max_degree' variables are None, one should be given")
    if not constr:
        return build_qp(n_vars, h, l)

    if SPP:
        A = rand_constr_matrix_SPP(n_vars, n_constr, max_degree = max_degree)
        b = np.ones(n_constr)
    else:
        if density is not None:
            A = rand_constr_matrix_density(n_vars, n_constr, density)
        elif max_degree is not None:
            A = rand_constr_matrix_degree(n_vars, n_constr, max_degree)
        else:
            raise ValueError("Both 'density' and 'max_degree' variables are None, one should be given")
        b = rand_constr_vector(n_constr)
    return build_qp(n_vars, h, l, A, b) 




def build_instance(normalize_method, n, n_cons, test_set, density, max_degree, seed = 42, to_file = True, SPP = False):
    '''
    Given the problem size (n_vars and n_cons), create a problem instance, normalize it, check feasibility and write to file
    normalize_method identifies the way in which the instance is generated and normalized. Specifically, BN, SN and NN
    refer to BruteNormalization, SamplingNormalization and NoNormalization, respectively.
    Return:
        p - the problem instance
    '''
    constr = True
    if n_cons == 0:
        constr = False

    # produce instance
    np.random.seed(seed)
    if normalize_method == "BN":
        qp = random_qp(n, constr, n_cons, density, max_degree)
        p = Problem(qp)
        H = p.get_obj_hamiltonian()
        norm_H = np.max(np.abs(H))
        H = H/norm_H
        p.normalize_problem(norm_H)
    elif normalize_method == "NN":
        qp = random_qp(n, constr, n_cons, density, max_degree, SPP = SPP)
        p = Problem(qp)
    elif normalize_method == "SN":
        qp = random_qp(n, constr, n_cons, density, max_degree, sampling_normalization = True)
        p = Problem(qp)
    else:
        raise ValueError(f"normalize_method {normalize_method} not among the implemented ones!")

    # test feasibility of instance
    # write
    file = open("test_feasibility.lp", "w")
    file.write(p.qp.export_as_lp_string())
    file.close()
    # read and try to solve
    m = ModelReader.read("test_feasibility.lp", ignore_names=True)
    p_test = Problem(from_docplex_mp(m))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        c_res = p_test.solve_exact()
    if c_res.status != OptimizationResultStatus.SUCCESS:
        print("Infeasible problem, trying another one")
        return build_instance(normalize_method, n, n_cons, test_set, density = density, max_degree = max_degree, seed = seed+1, to_file = to_file, SPP = SPP)
    if len(c_res.x) < n:
        print("Problem with less variables, trying another one")
        return build_instance(normalize_method, n, n_cons, test_set, density = density, max_degree = max_degree, seed = seed+1, to_file = to_file, SPP = SPP)
    osremove("./test_feasibility.lp")  # added line, still need to check if it works

    # write to LP file
    if to_file:
        filename = f"{test_set}/{n}/random{seed}_{n}_{n_cons}.lp"
        p.qp.write_to_lp_file(filename)
        print(f"File {filename} wrote")
    return p


def build_database(test_set, normalize_method, n_instances, n, n_cons = None, starting_seed = 42, density = None, max_degree = None, SPP = False):
    '''
    Write to file n_instances instances of problems to build a database
    '''
    if n_cons == None:
        n_cons = np.max([1, int(n/3)]) # divided by 3 for SPP, divided by 5 for linear

    seed = starting_seed
    for i in range(n_instances):
        build_instance(normalize_method, n, n_cons, test_set, density, max_degree, seed, SPP = SPP)
        seed += 100
    return


""" conn = 10
p = build_instance("NN", 6, 1, "ajo", None, conn, seed = 42, to_file = False)
print(p.qp.export_as_lp_string()) """

#print("\n\n\nChange hardcoded parameters in line 123 if you're changing degree or p_baseline\n\n\n")

# BUILD DATABASE
""" bvars = np.arange(6, 31)
degree = 5
n_samples = 1000
for nvars in bvars:
    build_database("../problems/SPP_p15", "NN", n_samples, nvars, max_degree = degree, SPP=True) """