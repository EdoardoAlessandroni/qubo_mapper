################ NN ################


NN_linear_deg10.txt
No Normalization of the Hamiltonian
n = range(4,26) + [30, 35, 40, 45, 50]
200 samples
Analysis of gap for Abiabatic evolution: (diag of mix hamiltonian) only in range(4,14).
Analysis of normalized gap (computing the hamiltonian) only in range(4,26).
Linear: number of terms in f is linear; Q matrix has density such that the every node has fixed max degree (10).
M_strategies = ["our_M", "qiskit_M", "optimal_M"]

NN_linear_deg5.txt
No Normalization of the Hamiltonian
n = range(4,26)
1000 samples
Analysis of normalized gap (computing the hamiltonian) only in range(4,26).
Linear: number of terms in f is linear; Q matrix has density such that the every node has fixed max degree (5).
M_strategies = ["our_M", "qiskit_M", "optimal_M", "babbush_M"]


################ SPP ################


SPP_p15.txt
Set Partitioning Problem instances
n = range(6,26)
200 samples
Analysis of normalized gap
p15: constraint matrix with density (prob any element is non zero) 0.15
M_strategies = ["our_M", "qiskit_M", "optimal_M"]

SPP_p25.txt
Set Partitioning Problem instances
n = range(6,26)
1000 samples
Analysis of normalized gap
p25: constraint matrix with density (prob any element is non zero) 0.25
M_strategies = ["our_M", "qiskit_M", "optimal_M", "babbush_M"]

SPP_p35.txt
Set Partitioning Problem instances
n = range(6,26)
200 samples
Analysis of normalized gap
p35: constraint matrix with density (prob any element is non zero) 0.35
M_strategies = ["our_M", "qiskit_M", "optimal_M"]


################ PO ################


PO_part3_mult2.txt
Portfolio Optimization, Markowitz model
n = range(6,25,3)
200 samples
pars: w=3, mult=1e2, ra=1
M_strategies = ["our_M", "qiskit_M", "optimal_M", "babbush_M"]

PO_part3_mult4.txt
Portfolio Optimization, Markowitz model
n = range(6,25,3)
1000 samples
pars: w=3, mult=1e4, ra=1
M_strategies = ["our_M", "qiskit_M", "optimal_M", "babbush_M"]

PO_norm_part3_mult4.txt
Portfolio Optimization, Markowitz model, (NORMalized formulation)
n = range(6,25,3)
1000 samples
pars: w=3, mult=1e4, ra=1
M_strategies = ["our_M", "qiskit_M"]

PO_norm_part3_mult2_test.txt
Portfolio Optimization, Markowitz model, (NORMalized formulation)
n = range(6,22,3)
100 samples
pars: w=3, mult=1e2, ra=1
M_strategies = ["our_M", "qiskit_M"]

PO_greedy
Portfolio Optimization, Markowitz model
n = range(6,25,3)
1000 samples
pars: w=3, mult=1e4, ra=1
M_strategies = ["heuristic_PO_M"]

PO_norm_big_greedy.txt
Portfolio Optimization, Markowitz model, (NORMalized formulation)
n = range(10,51,10) + [75] + range(100, 301, 50)
100 samples
pars: w=3, mult=1e4, ra=1
M_strategies = ["heuristic_PO_M", "qiskit_M"]


