
#### OPTIMALITY AND FEASIBILITY PERCENTAGES AND PLOT IT


def evaluate_optimality(data, opt_threshold = .5):
    n_bvars, n_samples, n_M_strats = np.shape(data.is_feasible)
    for var_idx in range(n_bvars):
        for M_idx in range(n_M_strats):
            for sample in range(n_samples):
                if data.is_feasible[var_idx, sample, M_idx]:
                    fc = data.fval[var_idx, sample, 0, M_idx]
                    fq = data.fval[var_idx, sample, 1, M_idx]
                    if np.abs(fc - fq) < opt_threshold: # optimum
                        data.is_optimum[var_idx, sample, M_idx] = True

def compute_errors(data, threshold_fc_zero = 1e-8):
    n_bvars, n_samples, n_M_strats = np.shape(data.is_feasible)
    is_feas_non_opt = np.logical_xor(data.is_feasible, data.is_optimum)
    for var_idx in range(n_bvars):
        for M_idx in range(n_M_strats):
            for sample in range(n_samples):
                if is_feas_non_opt[var_idx, sample, M_idx] == True:
                    fc = data.fval[var_idx, sample, 0, M_idx]
                    fq = data.fval[var_idx, sample, 1, M_idx]
                    if np.abs(fc) > threshold_fc_zero: 
                        ratio = np.abs( (fq -fc) / fc)
                    else:
                        ratio = None
                    data.relative_error[var_idx, sample, M_idx] = ratio
                    data.absolute_error[var_idx, sample, M_idx] = fq -fc

def statistics_correctness_n_feasability(data):
    '''
    Returns:
        percentage_optimal, percentage_feasible, (avg_relative_error, std_relative_error) [only among feasible non optimal solutions],
        (avg_n_violations, avg_max_violation) [only among non feasible solutions]
    '''
    n_bvars, n_samples, n_M_strats = np.shape(data.is_feasible)
    # percentages
    perc_opt = np.count_nonzero(data.is_optimum, axis=1) / n_samples
    perc_feas = np.count_nonzero(data.is_feasible, axis=1) / n_samples
    # relative errs
    is_feas_non_opt = np.logical_xor(data.is_feasible, data.is_optimum)
    avg_relative_error, std_relative_error = np.zeros((n_bvars, n_M_strats)), np.zeros((n_bvars, n_M_strats))
    for var_idx in range(n_bvars):
        for M_idx in range(n_M_strats):
            relative_errs = data.relative_error[var_idx, :, M_idx][is_feas_non_opt[var_idx, :, M_idx]]
            # discard the None for when fc ~ 0
            relative_errs = relative_errs[~np.isnan(relative_errs)]
            if len(relative_errs) > 0:
                avg_relative_error[var_idx, M_idx] = np.mean(relative_errs)
                std_relative_error[var_idx, M_idx] = np.std(relative_errs)
    # errors
    avg_error, std_error = np.zeros((n_bvars, n_M_strats)), np.zeros((n_bvars, n_M_strats))
    for var_idx in range(n_bvars):
        for M_idx in range(n_M_strats):
            absolute_err = data.absolute_error[var_idx, :, M_idx][is_feas_non_opt[var_idx, :, M_idx]]
            if len(absolute_err) > 0:
                avg_error[var_idx, M_idx] = np.mean(absolute_err)
                std_error[var_idx, M_idx] = np.std(absolute_err)
    # violations
    is_infeas = np.logical_not(data.is_feasible)
    avg_n_violations, avg_max_violation = np.zeros((n_bvars, n_M_strats)), np.zeros((n_bvars, n_M_strats))
    for var_idx in range(n_bvars):
        for M_idx in range(n_M_strats):
            n_viol = data.n_violations[var_idx, :, M_idx][is_infeas[var_idx, :, M_idx]]
            if len(n_viol) > 0:
                max_viol = data.max_violation[var_idx, :, M_idx][is_infeas[var_idx, :, M_idx]]
                avg_n_violations[var_idx, M_idx] = np.mean(n_viol)
                avg_max_violation[var_idx, M_idx] = np.mean(max_viol)
    return perc_opt, perc_feas, (avg_relative_error, std_relative_error), (avg_error, std_error), (avg_n_violations, avg_max_violation)

evaluate_optimality(data)
compute_errors(data)
perc_opt, perc_feas, (avg_relative_error, std_relative_error), (avg_absolute_error, std_absolute_error), (avg_n_violations, avg_max_violation) = statistics_correctness_n_feasability(data)



if np.any(data.is_optimum == 0):
    print("!!!!\tSomething was not optimally solved!")
    n_bvars, n_samples, n_Ms = np.shape(data.is_optimum)
    for i in range(n_bvars):
        print(f"Vars = {data.bvars[i]}\t solved {np.sum(data.is_optimum[i])} out of {n_samples*n_Ms}")
else:
    print("ALL was not optimally solved!")
    
    

fig = plt.figure(figsize=(3,2))
ax1 = fig.add_subplot(1,1,1)
ax1.set_title('Percentage optimality', fontsize = 14)
ax1.set_xlabel('# of binary variables', fontsize = 14)
ax1.set_ylabel('$\%$ optimal solution', fontsize = 14)
plt.ylim(0, 1.1)
for i in range(len(M_strategies)):
    if plot_M[i]:
        plt.plot(data.bvars, perc_opt[:,i], label=M_strategies[i])
plt.ylim(.75,1.02)
plt.grid()
ax1.legend(title = 'M choice', fontsize = 8, loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()

fig = plt.figure(figsize=(3,2))
ax1 = fig.add_subplot(1,1,1)
ax1.set_title('Percentage feasibility', fontsize = 14)
ax1.set_xlabel('# of binary variables', fontsize = 14)
ax1.set_ylabel('$\%$ feasible solution', fontsize = 14)
plt.ylim(.75,1.02)
for i in range(len(M_strategies)):
    if plot_M[i]:
        plt.plot(data.bvars, perc_feas[:,i], label=M_strategies[i])
plt.grid()
ax1.legend(title = 'M choice', fontsize = 8, loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()







### ERRORS IN SOLUTION FOUND: PLOTS

fig = plt.figure(figsize=(6,4))
ax1 = fig.add_subplot(1,1,1)
ax1.set_title('Relative error (among feas non-opt)', fontsize = 22)
ax1.set_xlabel('# of binary variables', fontsize = 22)
ax1.set_ylabel('$|f_q-f_c / f_c|$', fontsize = 22)
is_feas_non_opt = np.logical_xor(data.is_feasible, data.is_optimum)
for i in range(len(M_strategies)):
    if plot_M[i]:
        plt.plot(data.bvars, avg_relative_error[:,i], label=M_strategies[i], color = colors[i])
        plt.fill_between(data.bvars, avg_relative_error[:,i] - std_relative_error[:,i], avg_relative_error[:,i] + std_relative_error[:,i], alpha=.2)
#plt.ylim(-1,5)
plt.grid()
ax1.legend(title = 'M choice', fontsize = 14, loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()

fig = plt.figure(figsize=(6,4))
ax1 = fig.add_subplot(1,1,1)
ax1.set_title('Absolute error (among feas non-opt)', fontsize = 22)
ax1.set_xlabel('# of binary variables', fontsize = 22)
ax1.set_ylabel('$f_q-f_c$', fontsize = 22)
for i in range(len(M_strategies)):
    if plot_M[i]:
        plt.plot(data.bvars, avg_absolute_error[:,i], label=M_strategies[i], color = colors[i])
        plt.fill_between(data.bvars, avg_absolute_error[:,i] - std_absolute_error[:,i], np.abs(avg_absolute_error[:,i]) + std_absolute_error[:,i], alpha=.2)
#plt.ylim(0,2)
plt.grid()
ax1.legend(title = 'M choice', fontsize = 14, loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()

fig = plt.figure(figsize=(6,4))
ax1 = fig.add_subplot(1,1,1)
ax1.set_title('Average number violations (among non feasible)', fontsize = 22)
ax1.set_xlabel('# of binary variables', fontsize = 22)
for i in range(len(M_strategies)):
    if plot_M[i]:
        plt.plot(data.bvars, avg_n_violations[:,i], label=M_strategies[i])
plt.grid()
ax1.legend(title = 'M choice', fontsize = 14, loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()

fig = plt.figure(figsize=(6,4))
ax1 = fig.add_subplot(1,1,1)
ax1.set_title('Average maximum violation (among non feasible)', fontsize = 22)
ax1.set_xlabel('# of binary variables', fontsize = 22)
for i in range(len(M_strategies)):
    if plot_M[i]:
        plt.plot(data.bvars, avg_max_violation[:,i], label=M_strategies[i])
plt.grid()
ax1.legend(title = 'M choice', fontsize = 14, loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()








### CHERRYPICKING BIG GAP INSTANCES


import numpy as np
import matplotlib.pyplot as plt
from ds import Datas

var_idx = -2
gaps_o = data.gap_norm[var_idx,:,0]
idx = np.argsort(gaps_o)

fig = plt.figure(figsize=(6,4))
ax1 = fig.add_subplot(1,1,1)
plt.plot(np.arange(len(gaps_o)), gaps_o[idx], ".", color = "C0")
plt.plot(np.arange(len(gaps_q)), gaps_q[idx], ".", color = "C1")
plt.grid()
plt.yscale("log")
plt.show()

# cherrypicking 
size_cherries = 25
ind = np.argpartition(gaps_o, -size_cherries)[-size_cherries:]

fig = plt.figure(figsize=(6,4))
ax1 = fig.add_subplot(1,1,1)
plt.plot(np.arange(len(ind)), gaps_o[ind], ".", color = "C0")
plt.plot(np.arange(len(ind)), gaps_q[ind], ".", color = "C1")
plt.grid()
plt.yscale("log")
plt.show()

#f = open("../easy_toys_adiabevol/PO_sp500_part3_ra10_mult2/easy_inst21.txt", "w")
#for i in ind:
#   f.write(data.filenames[var_idx,i]+"\n")
#f.close()







### TESTING HYPERPLANE ROUNDING ON AN INSTANCE

import cvxpy as cp
from qiskit_optimization.translators import from_docplex_mp
from docplex.mp.model_reader import ModelReader
from ds import Problem
from qiskit_optimization import QuadraticProgram
from scipy.linalg import cholesky, eigvalsh

def hyperplane_rounding(Y, add_on = 1e-7):
    n, _ = np.shape(Y)
    #print(f"Y:\t{Y}")
    Y = Y + add_on*np.eye(n)
    #eig = eigvalsh(Y)
    #print(eig)

    L = cholesky(Y)
    n_samp = 100
    g = np.random.randn(n_samp, n)
    x = (np.dot(g, L) > 0).astype(int)
    x = np.rint( np.mean(x, axis = 0) ).astype(int)

    if x[0] == 1:
        return x[1:]
    return (np.ones(n-1) - x[1:] ).astype(int)

def solve_SDP(p):
    n = p.n_vars
    Q = p.obj_quadratic
    L = p.obj_linear
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
    sol = prob.variables()[0].value
    return sol, prob.value

def solve(p):
    n = p.n_vars
    Q = p.obj_quadratic
    L = p.obj_linear
    qp = QuadraticProgram()
    for i in range(n):
        qp.binary_var()
    qp.minimize(linear = L, quadratic = Q)
    p_unc = Problem(qp)
    result = p_unc.solve_exact()
    return result.x, result.fval

filename = "../toys/NN_linear_deg5/18/random10842_18_3.lp"
m = ModelReader.read(filename, ignore_names=True)
qp = from_docplex_mp(m)
p = Problem(qp)

Y, fval_sdp = solve_SDP(p)
x, fval = solve(p)
x = x.astype(int)

x_hr = hyperplane_rounding(Y)
print(f"fval\nsdp:\t{fval_sdp}\treal:\t{fval}\tapprox:\t{p.qp.objective.evaluate(x_hr)}")
print(f"points\nappr\t{x_hr}\nglobal\t{x}")









### TESTING HYPERPLANE ROUNDING ON A DATASET

import cvxpy as cp
from qiskit_optimization.translators import from_docplex_mp
from docplex.mp.model_reader import ModelReader
from ds import Problem
from qiskit_optimization import QuadraticProgram
from scipy.linalg import cholesky, eigvalsh
from os import listdir
from time import time

def run_instance_HR(filename):
    '''
    Read LP file to get problem instance, solve it both in a classic and quantum way(s) and compute the gaps
    Return:
        p - the problem instance
        xs - the results got with the M strategies
    '''
    m = ModelReader.read(filename, ignore_names=True)
    qp = from_docplex_mp(m)
    p = Problem(qp)
    p.qp.name = filename
    
    n = p.n_vars
    Q = p.obj_quadratic
    L = p.obj_linear
    qp = QuadraticProgram()
    for i in range(n):
        qp.binary_var()
    qp.minimize(linear = L, quadratic = Q)
    p_unc = Problem(qp)
    c_res = p_unc.solve_exact()
    x_global = np.rint(c_res.x).astype(int)
    
    _, Y = p.solve_unconstrained("SDP")
    # Hyperplane Rounding
    add_on = 1e-7
    n_samp = 100

    n, _ = np.shape(Y)
    Y = Y + add_on*np.eye(n)
    #eig = eigvalsh(Y)
    #print(eig)
    L = cholesky(Y)
    g = np.random.randn(n_samp, n)
    x = (np.dot(g, L) > 0).astype(int)
    x = np.rint( np.mean(x, axis = 0) ).astype(int)

    if x[0] == 1:
        x_hr =  x[1:]
    else:
        x_hr = (np.ones(n-1) - x[1:] ).astype(int)

    n_flips = np.sum(np.abs(x_global - x_hr)).astype(int)
    m = len(p.constraints)
    violations = np.ndarray((2, m))
    for i in range(m):
        cons = p.constraints[i]
        violations[0,i] = (cons.evaluate(x_global) - cons.rhs)**2
        violations[1,i] = (cons.evaluate(x_hr) - cons.rhs)**2

    return n_flips, violations


def run_test(test_set, bvars, n_samples):
    '''
    Run simulation of problems (read from files) for different number of qubits, M-choice strategies, and samples and return data acquired
    '''
    if test_set[8:10] == "SP":
        m_max = (bvars[-1]/3).astype(int)
    elif test_set[8:10] == "NN":
        m_max = (bvars[-1]/5).astype(int)
    else:
        raise ValueError("What kind of test dataset is even that?!")
    print(m_max)

    flips = np.ndarray((len(bvars), n_samples))
    violat = np.ndarray((len(bvars), n_samples, 2, m_max))

    for i in range(len(bvars)):
        n_qubs = bvars[i]
        m = (n_qubs/5).astype(int)
        print("\n" + str(n_qubs))
        folder = test_set+"/"+str(n_qubs)+"/"
        files = sorted(listdir(folder))
        if len(files) < n_samples:
            raise ValueError(f"Folder {folder} contains only {len(files)} instances, {n_samples} were requested")
        
        tic = time()
        for sample in range(n_samples):
            filename = folder + files[sample]
            print(sample, end = ", ")
            n_flips, violations = run_instance_HR(filename)
            flips[i, sample] = n_flips
            violat[i, sample] = np.pad(violations, [(0,0),(0, m_max - m)])
        tac = time()
        print(f"It took {tac - tic} sec")
    return flips, violat


bvars = np.arange(6, 23, 2)
n_samples = 200
test_set = "../toys/NN_linear_deg5"
flips, ks = run_test(test_set, bvars, n_samples)









### PLOTTING RESULTS OF HYPERPLANE ROUNDING TEST ON A DATASET

avg, std = np.mean(flips, axis = 1), np.std(flips, axis = 1)

fig = plt.figure(figsize=(6,4))
ax1 = fig.add_subplot(1,1,1)
ax1.set_title('Percentage of errors', fontsize = 22)
ax1.set_xlabel('n', fontsize = 22)
ax1.set_ylabel('wrong_bits / n', fontsize = 22)
plt.plot(bvars, avg/bvars)
plt.fill_between(bvars, avg/bvars - std/bvars, avg/bvars + std/bvars, alpha=.2)
plt.grid()
plt.show()

fig = plt.figure(figsize=(6,4))
ax1 = fig.add_subplot(1,1,1)
ax1.set_title('Errors', fontsize = 22)
ax1.set_xlabel('n', fontsize = 22)
ax1.set_ylabel('wrong_bits', fontsize = 22)
plt.plot(bvars, avg)
plt.fill_between(bvars, avg - std, avg + std, alpha=.2)
plt.grid()
plt.show()


k_gl = np.sum(ks[:, :, 0], axis = -1)
k_hr = np.sum(ks[:, :, 1], axis = -1)

avg = np.ndarray((len(bvars)))
std = np.ndarray((len(bvars)))
for i in range(len(bvars)):
    idx = 0
    l = []
    for s in range(n_samples):
        if k_gl[i, s] != 0:
            idx += 1
            l.append(np.abs(k_gl[i, s] - k_hr[i, s])/k_gl[i, s])
    avg[i] = np.mean(l)
    std[i] = np.std(l)


fig = plt.figure(figsize=(6,4))
ax1 = fig.add_subplot(1,1,1)
ax1.set_title('Relative error on violation number', fontsize = 22)
ax1.set_xlabel('n', fontsize = 22)
ax1.set_ylabel('|k - k_hr| / k', fontsize = 22)
plt.plot(bvars, avg)
#plt.fill_between(bvars, avg - std, avg + std, alpha=.2)
plt.ylim(0, 2)
plt.grid()
plt.show()