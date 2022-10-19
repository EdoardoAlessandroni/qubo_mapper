BN.txt
Results for n in range(2,24,2) and bruteforce normalization and gap analysis
M_strategies = ["qiskit_M", "our_M", "naive_M", "our_M*0.5"]

SN.txt
Results for n in range(2,21,2) + [24, 30, 35, 40] with normalization by sampling objective function coefficients quadratically wrt inverse of number of bvars and no gap analysis
M_strategies = ["qiskit_M", "our_M", "naive_M", "our_M*0.5"]

NN.txt
== data_NN_gap.txt + data_NN_nogap.txt
Results for n in range(2,25,2) + [30, 35] with no normalization of the problem and with gap analysis only in range(2,25,2)
M_strategies = ["qiskit_M", "our_M", "naive_M"]

NN_50.txt
== data_NN_50_allgaps.txt + data_NN_50_partialgaps.txt + data_NN_50_nogaps.txt
Results for n in range(4,25,2) + [30, 35] with no normalization of the problem, with gap analysis for Abiabatic evolution (requires diagonalization of mixture hamiltonian) only in range(4,13,2) and with gap analysis for QITE evolution (requires computing the hamiltonian) only in range(4,25,2)
50 stands for the fact that Q matrix has 0.5 density
M_strategies = ["our_M", "qiskit_M", "optimal_M"]

NN_linear_deg10.txt
== data_NN_linear_fit_allgaps.txt + data_NN_linear_fit_partialgaps.txt + data_NN_linear_fit_nogaps.txt
Results for n in range(4,26) + [30, 35, 40, 45, 50] with no normalization of the problem, with 200 samples (rather than 100 as usual).
With gap analysis for Abiabatic evolution (requires diagonalization of mixture hamiltonian) only in range(4,14) and with gap analysis for QITE evolution (requires computing the hamiltonian) only in range(4,26).
Linear stands for the fact that Q matrix has density such that the every node has fixed max degree (10), therefore the number of terms in f is linear and the density of Q is inversely proportional to n.
M_strategies = ["our_M", "qiskit_M", "optimal_M"]

NN_linear_deg5.txt
Results for n in range(4,26) with no normalization of the problem, with 200 samples (rather than 100 as usual).
With gap analysis for QITE evolution (requires computing the hamiltonian) in full range(4,26).
Linear stands for the fact that Q matrix has density such that the every node has fixed max degree (5), therefore the number of terms in f is linear and the density of Q is inversely proportional to n.
M_strategies = ["our_M", "qiskit_M", "optimal_M"]

SPP_deg5.txt
Results for n in range(6,26) with no normalization of the problem, with 200 samples.
With gap analysis for QITE evolution (requires computing the hamiltonian) in full range(6,26).
SPP stands for the fact that Set Partitioning problems are analyzed (from SPP_deg5). A matrix has density such that the every node has fixed max mean degree (5).
M_strategies = ["our_M", "qiskit_M"]

SPP_deg10.txt
Results for n in range(6,26) with no normalization of the problem, with 200 samples.
With gap analysis for QITE evolution (requires computing the hamiltonian) in full range(6,26).
SPP stands for the fact that Set Partitioning problems are analyzed (from SPP_deg10). A matrix has density such that the every node has fixed max mean degree (5).
M_strategies = ["our_M", "qiskit_M"]