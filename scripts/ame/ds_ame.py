import numpy as np
import cvxpy as cp
from qiskit_optimization.converters import LinearEqualityToPenalty
from qiskit_optimization import QuadraticProgram
from qiskit_optimization.algorithms import CplexOptimizer
import cplex
import time


def initial_hamiltonian_adiabatic(n):
    H = -kron_x_i(0,n)
    for i in range(1, n):
        H -= kron_x_i(i, n)
    return H

def kron_x_i(i, n):
    x = np.array([[0,1],[1,0]])
    id1 = np.eye(int(2**(i)))
    id2 = np.eye(int(2**(n-i-1)))
    return np.kron( np.kron(id1, x), id2 )


class Datas():
    def __init__(self, bvars, n_samples):
        self.bvars = bvars
        n_bvars = len(bvars) # bvars is a list containing the number of binary variables we wish to investigate, not necessarily consecutive
        self.is_optimum = np.zeros((n_bvars, n_samples), dtype = bool)
        #self.is_feasible = np.zeros((n_bvars, n_samples), dtype = bool)
        #self.time = np.ndarray((n_bvars, n_samples))
        self.fval_classic = np.ndarray((n_bvars, n_samples), dtype = int)
        self.max_iter = np.ndarray((n_bvars, n_samples), dtype = int)

        # in following dictionaries, we frist recover the data by accessing with the key n_bvars_n_sample
        #[n_bvars, n_samples, max_iter_tot]
        self.gaps = {} # gap referring to    H_ob + M H_c shifted and squeezed
        self.fvals = {}
        self.Ms = {}
        self.violation_nums = {}
        

# class Problem to tie together the QuadraticProgram from qiskit, its Ising formulation, the Ising formualtion of the constraints and their respective (obj and constraints) Hermitian matrix representation 2^n x 2^n.
# It also deals with solving the problem with different M strategies
class Problem():
    def __init__(self, quadratic_problem):

        self.n_vars = quadratic_problem.get_num_binary_vars()
        self.qp = quadratic_problem

        # objective
        self.obj_linear = quadratic_problem.objective.linear.to_array()
        self.obj_quadratic = quadratic_problem.objective.quadratic.to_array()

        # constraints
        m = quadratic_problem.get_num_linear_constraints()
        if m == 0:
            self.constrained = False
        else:
            self.constrained = True
            self.constraints = []
            for i in range(m):
                self.constraints.append( quadratic_problem.get_linear_constraint(i) )
    

    def get_constr_matrices(self):
        '''
        Get the A matrix and the b vector of the constraints from the QuadraticProgram attributes
        '''
        m = self.qp.get_num_linear_constraints()
        n = self.n_vars
        A = np.zeros((m,n))
        b = np.zeros((m))       
        
        for i in range(m):
            b[i] = self.constraints[i].rhs
            coefs = self.constraints[i].linear.coefficients
            keys = list(coefs.keys())
            for col in keys:
                A[i, col[1]] = coefs[col]
        return A, b


    def to_ising(self):
        '''
        Get the J (triangular) and h matrices of the Ising formulation from the Q and L matrices of the QUBO formulation
        '''
        # first we move the terms of the form q_ii x_i x_i to the linear vector as x_i^2 = x_1
        L = self.obj_linear + np.diag(self.obj_quadratic)
        Q = self.obj_quadratic
        n = self.n_vars
        # this loop is not needed since we don't access Q_ii anymore
        #for i in range(n):
        #    Q[i,i] = 0
        
        # then we map it to the J and h matrix of the Ising formulation
        h = np.ndarray(n)
        J = np.zeros((n, n))
        const_term = np.sum(L)/2

        for i in range(n):
            h[i] = L[i]/2 + (np.sum(Q[i, i+1:]) + np.sum(Q[:i, i]))/4
            const_term += np.sum(Q[i, i+1:])/4
            for j in range(i+1, n):
                J[i,j] = Q[i,j]/4
        return J, h


    def constraints_to_qubo_form(self):
        '''
        Get the Q matrix and L vector of the QUBO formualtion from the QuadraticProgram
        '''
        if self.constrained == False:
            raise ValueError("The problem is unconstrained, no constraints to map to QUBO")
        A, b = self.get_constr_matrices()

        Q = A.T.dot(A)
        # make Q triangular
        Q = np.triu(Q + Q.T - np.diag(np.diag(Q)))
        L = -2*A.T.dot(b)
        const_term = b.T.dot(b) # not returned but we compute it if we need it in the future
        return Q, L, const_term


    def constraints_to_ising_form(self):
        '''
        Get the J (triangular) matrix and h vector of the ising formualtion from the QuadraticProgram
        '''
        Q, L, old_const_term = self.constraints_to_qubo_form()
        n = self.n_vars

        # move terms from diagonal
        L += np.diag(Q)
        # this loop is not needed since we don't access Q_ii anymore
        #for i in range(self.n_vars):
        #    Q[i,i] = 0
    
        # map to ising formulation
        h = np.ndarray(n)
        J = np.zeros((n, n))
        const_term = np.sum(L)/2 + old_const_term

        for i in range(n):
            h[i] = L[i]/2 + (np.sum(Q[i, i+1:]) + np.sum(Q[:i, i]))/4
            const_term += np.sum(Q[i, i+1:])/4
            for j in range(i+1, n):
                J[i,j] = Q[i,j]/4
        return J, h, const_term

    
    def from_ising_to_hamiltonian(self, J, h, const = 0):
        '''
        Get the 2^nX2^n Hamiltonian corresponding to a given Ising formulation of the nXn matrix J and vector h
        '''
        n = self.n_vars
        H = np.zeros((2**n))
        # coupling terms
        for i in range(n):
            for j in range(i+1, n):
                if np.abs(J[i,j]) >= 1e-12:
                    H += J[i,j]*self.__kron_ij__(i,j)
        # field terms
        for i in range(n):
            if np.abs(h[i]) >= 1e-12:
                H += h[i]*self.__kron_i__(i)
        return H + const

    
    def __kron_ij__(self, i, j):
        '''
        Compute tensor product 
        id X .. X sigma_i X .. X sigma_j X .. X id
        '''
        n = self.n_vars
        res = np.ones((1))
        for k in range(n-1, -1, -1):
            if k != i and k != j:
                res = np.concatenate([res, res])
            else:
                res = np.concatenate([res, -res])
        return res
        

    def __kron_i__(self, i):
        '''
        Compute the diagonal of the tensor product 
        id X .. X sigma_i X .. X id
        '''
        
        n = self.n_vars
        res = np.ones((1))
        for k in range(n-1, -1, -1):
            if k != i:
                res = np.concatenate([res, res])
            else:
                res = np.concatenate([res, -res])
        return res


    def get_obj_hamiltonian(self):
        J, h = self.to_ising()
        return self.from_ising_to_hamiltonian(J, h)
    
    
    def get_constraint_hamiltonian(self):
        J, h, const = self.constraints_to_ising_form()
        return self.from_ising_to_hamiltonian(J, h, const)
    

    def normalize_problem(self, norm_H):
        # normalize objective function
        self.obj_quadratic = self.obj_quadratic / norm_H
        self.obj_linear = self.obj_linear / norm_H
        self.qp.minimize(linear = self.obj_linear, quadratic = self.obj_quadratic)
        return


    def get_gap_objective(self, evs_H):
        evs = np.unique(evs_H)
        #evs = np.partition(evs, kth=1)[:2]
        return evs[1] - evs[0]
    

    def get_gap_total(self, H, Hc, M):
        evs = H + M*Hc
        evs = np.unique(evs) # already sorted
        #evs = np.partition(evs, kth=1)[:2]
        return evs[1] - evs[0]
    
    def solve_exact(self):
        result = CplexOptimizer(disp = False).solve(self.qp)
        return result

    def solve_quantum(self, M_strategy):
        # choose M
        if M_strategy == "qiskit_M":
            converter = LinearEqualityToPenalty()
        elif M_strategy == "our_M":
            M = self.our_M()
            converter = LinearEqualityToPenalty(penalty = M)
        elif M_strategy == "naive_M":
            M = self.naive_M()
            converter = LinearEqualityToPenalty(penalty = M)
        elif M_strategy == "optimal_M":
            M = self.optimal_M()
            converter = LinearEqualityToPenalty(penalty = M)
        elif M_strategy[:6] == "our_M*":
            multiplier = float(M_strategy[6:])
            M = self.our_M()*multiplier
            converter = LinearEqualityToPenalty(penalty = M)
        else:
            raise ValueError(f"M_strategy {M_strategy} not known")
        qubo = converter.convert(self.qp)
        M = converter.penalty
        
        tic = time.time()
        result = CplexOptimizer(disp = False).solve(qubo)
        tac = time.time()
        return result, M, tac - tic


    def naive_M(self): # adds all positive coefficients in Q matrix and L vector of the objective function
        L = self.obj_linear
        Q = self.obj_quadratic
        return np.sum(L*(L>0)) + np.sum(Q*(Q>0)) - np.sum(L*(L<0)) - np.sum(Q*(Q<0))

    def optimal_M(self):
        # find and evaluate feasible solution
        f = self.solve_exact().fval
        # evaluate unconstrained and continuous problem
        f_unc = self.solve_unconstrained(how = "exact")
        return f - f_unc + .5

    def our_M(self):
        # find and evaluate feasible solution
        f_feas = self.get_feasible_sol_objective()
        # evaluate unconstrained and continuous problem
        f_unc = self.solve_unconstrained(how = "SDP")
        return f_feas - f_unc + .5


    def solve_unconstrained(self, how = "SDP"):
        '''
        Solve the unconstrained relaxation of the problem, either with SDP (and MOSEK solver) or by getting exact solution (with CPLEX solver)
        '''
        n = self.n_vars
        Q = self.obj_quadratic
        L = self.obj_linear
        if how == "SDP":
            Q_tilde = np.ndarray((n+1, n+1))
            Q_tilde[1:, 1:] = Q
            Q_tilde[0, 0] = 0
            Q_tilde[0, 1:] = .5*L.T
            Q_tilde[1:, 0] = .5*L
            X = cp.Variable((n+1, n+1), symmetric=True)
            constraints = [X >> 0]
            constraints += [X[i,i] == X[0,i] for i in range(1, n+1) ]
            constraints += [X[i,j] <= 1 for i in range(n+1) for j in range(i, n+1) ]
            constraints += [X[i,j] >= 0 for i in range(n+1) for j in range(i, n+1) ]
            prob = cp.Problem(cp.Minimize(cp.trace(Q_tilde @ X)), constraints)
            prob.solve(solver = "MOSEK")
            return prob.value
        elif how == "exact":
            # build unconstrained problem
            qp = QuadraticProgram()
            for i in range(n):
                qp.binary_var()
            qp.minimize(linear = L, quadratic = Q)

            p_unc = Problem(qp)
            result = p_unc.solve_exact()
            return result.fval
        else:
            raise ValueError("How unconstrained relaxation should be solved is not among the implemented possibilities")

    
    def get_feasible_sol_objective(self):
        model = cplex.Cplex(self.qp.name)
        model.parameters.timelimit.set(10)
        model.set_results_stream(None)
        model.set_log_stream(None)
        model.solve()
        if not model.solution.is_primal_feasible():
            print("Using feasibility pump and mip emphasis for the first time!")
            model = cplex.Cplex(self.qp.name)
            model.parameters.mip.strategy.fpheur.set(1)
            model.parameters.emphasis.mip.set(1)
            model.parameters.timelimit.set(10)
            model.set_results_stream(None)
            model.set_log_stream(None)
            model.solve()
            if not model.solution.is_primal_feasible():
                print("Even by using feasibility pump and mip emphasis a feasible solution was not found")
        return model.solution.get_objective_value()


    def write_to_lp_file(self, filename = None):
        lp_string = self.qp.export_as_lp_string()
        if filename is None:
            filename = self.qp.name
        file = open(filename, "w")
        file.write(lp_string)
        file.close()
        return