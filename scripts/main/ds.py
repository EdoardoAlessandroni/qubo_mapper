import numpy as np
import cvxpy as cp
from qiskit_optimization.converters import LinearEqualityToPenalty
from qiskit_optimization import QuadraticProgram
from qiskit_optimization.algorithms import CplexOptimizer
import cplex
import time
from copy import deepcopy

# Classes built here: Datas, Problem


# get relevant data over many instances and bumber of variables with class Datas
class Datas():
    def __init__(self, bvars, n_samples, M_strategies):
        self.bvars = bvars
        self.M_strategies = M_strategies
        n_M_strategies = len(M_strategies)
        n_bvars = len(bvars) # bvars is a list containing the number of binary variables we wish to investigate, not necessarily consecutive
        self.filenames = np.ndarray((n_bvars, n_samples), dtype=object)
        # gaps are multiplied by -100 to spot problems: if we dont compute gaps but we plot it it's clear that we shouldn't look at it since that's negative
        self.gap_norm = -100*np.ones((n_bvars, n_samples, n_M_strategies)) # gap referring to    [H_ob + M H_c] / norm(H_ob + M H_c)  useful for QITE approach
        self.fval = np.ndarray((n_bvars, n_samples, 2, n_M_strategies), dtype = int) # 3rd index: 0->classical solution  1->quantum solution
        self.M = np.ndarray((n_bvars, n_samples, n_M_strategies))
        self.is_optimum = np.zeros((n_bvars, n_samples, n_M_strategies), dtype = bool)
        self.is_feasible = np.zeros((n_bvars, n_samples, n_M_strategies), dtype = bool)
        self.relative_error = np.zeros((n_bvars, n_samples, n_M_strategies))
        self.absolute_error = np.zeros((n_bvars, n_samples, n_M_strategies))
        self.n_violations = np.zeros((n_bvars, n_samples, n_M_strategies), dtype = int)
        self.max_violation = np.zeros((n_bvars, n_samples, n_M_strategies), dtype = int)


# class Problem to tie together the QuadraticProgram from qiskit, its Ising formulation, the Ising formualtion of the constraints and their respective (obj and constraints) Hamiltonian matrix representation 2^n x 2^n.
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
        Q = self.obj_quadratic - np.diag(np.diag(self.obj_quadratic))
        const = self.qp.objective.constant
        n = self.n_vars

        # then we map it to the J and h matrix of the Ising formulation
        h = np.ndarray(n)
        J = np.zeros((n, n))
        const_term = np.sum(L)/2 + const

        for i in range(n):
            h[i] = L[i]/2 + (np.sum(Q[i, i+1:]) + np.sum(Q[:i, i]))/4
            const_term += np.sum(Q[i, i+1:])/4
            for j in range(i+1, n):
                J[i,j] = Q[i,j]/4

        return J, h, const_term


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
        L = -2*b.T.dot(A) # it's equivalent, but before it was -2*A.T.dot(b)
        const_term = b.T.dot(b) # not returned but we compute it if we need it in the future
        # move diagonal terms of Q into L
        L += np.diag(Q)
        np.fill_diagonal(Q, 0)
        return Q, L, const_term


    def constraints_to_ising_form(self):
        '''
        Get the J (triangular) matrix and h vector of the ising formualtion from the QuadraticProgram
        '''
        Q, L, old_const_term = self.constraints_to_qubo_form()
        n = self.n_vars

        # move terms from diagonal
        L += np.diag(Q)            
        Q -= np.diag(np.diag(Q))   

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
        '''
        Get Ising Hamiltonian of the objective function 
        '''
        J, h, const = self.to_ising()
        return self.from_ising_to_hamiltonian(J, h, const)
    
    
    def get_constraint_hamiltonian(self):
        '''
        Get Ising Hamiltonian of the constraint penalization part 
        '''
        J, h, const = self.constraints_to_ising_form()
        return self.from_ising_to_hamiltonian(J, h, const)


    def get_gap_objective(self, evs):
        '''
        Get (non-normalized) gap of the input set of eigenvalues
        '''
        evs = np.unique(evs.round(decimals=14))
        return evs[1] - evs[0]
    

    def get_gap_total(self, H, Hc, M):
        '''
        Get (non-normalized) gap of the input Ising form in the penalized qubo formulation 
        '''
        evs = H + M*Hc
        evs = np.unique(evs.round(decimals=14))
        return evs[1] - evs[0]
    

    def solve_exact(self):
        '''
        Exactly solve the problem with CPLEX
        '''
        result = CplexOptimizer(disp = False).solve(self.qp)
        return result


    def solve_quantum(self, M_strategy):
        '''
        Solve the problem, mapping to the QUBO formulation with the specified M strategy.
        Returns the 'result' of CPLEX, together with M computed
        '''
        # choose M
        if M_strategy == "qiskit_M":
            converter = LinearEqualityToPenalty()
        elif M_strategy == "our_M":
            M = self.our_M()
            converter = LinearEqualityToPenalty(penalty = M)
        elif M_strategy == "optimal_M":
            M = self.optimal_M()
            converter = LinearEqualityToPenalty(penalty = M)
        elif M_strategy == "babbush_M":
            M = self.babbush_M()
            converter = LinearEqualityToPenalty(penalty = M)
        elif M_strategy == "heuristic_PO_M":
            M = self.heuristic_PO_M()
            converter = LinearEqualityToPenalty(penalty = M)
        else:
            raise ValueError(f"M_strategy {M_strategy} not known")
        
        qubo = converter.convert(self.qp)
        M = converter.penalty
        result = CplexOptimizer(disp = False).solve(qubo)
        return result, M
    

    def integer_mapper_PO(self, bin_per_int, mu, sigma, sigma_triangular = True):
        '''
        Maps the binary-suitable instance parameters (mu and sigma) to the integer-suitable form
        '''
        mu2 = mu[bin_per_int-1 :: bin_per_int]
        sigma2 = sigma[bin_per_int-1 :: bin_per_int, bin_per_int-1 :: bin_per_int]
        if sigma_triangular:
            sigma2 = (sigma2 + sigma2.T)/2
        return mu2, sigma2


    def greedy_heuristic(self, N, w, mu, sigma):
        '''
        Implements greedy algorithm to return nearly-optimal solution and its obj function value, given the instance parameters
        '''
        X = np.zeros((N), dtype = int)
        obj = np.empty((N))
        for j in range(2**w-1):
            for k in range(N):
                # compute new objective value, if we would add one unit (the j-th) of asset k
                new_X = deepcopy(X)
                new_X[k] += 1
                obj[k] = -np.dot(mu, new_X) + np.dot(new_X, sigma@new_X)
            winner = np.argmin(obj)
            X[winner] += 1
        f_X = -np.dot(mu, X) + np.dot(X, sigma@X)
        
        assert np.sum(X) == 2**w-1
        return X, f_X


    def heuristic_PO_M(self):
        '''
        Computes the greedy M: lower bound with SDP and upper bound with greedy PortOpt strategy (function: greedy_heuristic)
        '''
        # compute M by greedy heuristic
        n = self.n_vars # binary variables
        w = self.infer_partition_number()
        N = int(np.rint(n/w)) # integer variables

        sigma = self.obj_quadratic
        mu = -self.obj_linear
        mu, sigma = self.integer_mapper_PO(w, mu, sigma)
        
        # find and evaluate feasible solution
        X, f_feas = self.greedy_heuristic(N, w, mu, sigma)

        # evaluate unconstrained and continuous problem
        f_unc = self.solve_unconstrained(how = "SDP")
        return f_feas - f_unc + .5
    

    def infer_partition_number(self):
        '''
        Returns the partition number w of the instance, from its filename (contains "part[w]")
        '''
        filename = self.qp.name
        for idx in range(0, len(filename)):
            if filename[idx] == "p" and filename[idx+1 : idx+4] == "art":
                w = int(filename[idx+4])
                break

        if idx >= len(filename) -1:
            raise ValueError("Partition number can't be extracted from filename")
        return w


    def babbush_M(self):
        # compute M by max ( sum of all positive coeff,  sum of all negative coeff ) in obj funct
        Q = self.obj_quadratic
        L = self.obj_linear
        n = self.n_vars
        sum_neg = np.sum([Q[i,j] for i in range(n) for j in range(n) if Q[i,j] < 0]) + np.sum([L[i] for i in range(n) if L[i] < 0])
        sum_pos = np.sum([Q[i,j] for i in range(n) for j in range(n) if Q[i,j] > 0]) + np.sum([L[i] for i in range(n) if L[i] > 0])
        sum = np.max([sum_pos, -sum_neg])
        return 1 + sum


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
        #print(f"SDP relaxation took {tac - tic} seconds")
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
            #constraints += [X[0,0] == 1]      ###### alternative way of imposing boundness, rather than  0 <= X_ij <= 1  for all i,j (2 previous lines)
            prob = cp.Problem(cp.Minimize(cp.trace(Q_tilde @ X)), constraints)
            prob.solve(solver = "MOSEK")
            #X = prob.variables()[0].value
            return prob.value #, X
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
        '''
        Geta feasible solution, first running on CPLEX for a max of 10 secs, then, if the solution found is not feasible, another run starts with a focus on feasibility
        '''
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
        '''
        Export the Problem as a lp string and writes on file
        '''
        lp_string = self.qp.export_as_lp_string()
        if filename is None:
            filename = self.qp.name
        file = open(filename, "w")
        file.write(lp_string)
        file.close()
        return