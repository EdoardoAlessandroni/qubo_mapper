The folder contains QUBO formulation using different Ms (our method, qiskit method), alongside with constrained 
formulations and solutions (both x and f(x)), for "easy instances", meaning they have a big relative gap, and [6, 12, 18, 24] qubits.
10 instances for each dataset are present.

DIRECTORY STRUCTURE:

1st layer:
	NN_linear_deg5              --> unstructured dataset
	SPP_p25                     --> Set Partitioning Problem dataset
	PO_norm_part3_mult4   --> Portfolio Optimization dataset
	
2nd layer:
	n_qubits --> size of problem instance

3rd layer:
	inst_names.txt --> filenames of 25 big-gap instances
	constrained --> original constrained quadratic problem
	solution --> solution point (x) and value (f(x)) of the constrained problem
	ourM --> QUBO reformualtion using our method
	qiskM --> QUBO reformualtion using qiskit method

4th layer (only for qiskM and ourM 3rd layer directories):
	QUBO --> qubo reformulation
	gap --> normalized spectral gap
	maxener --> maximum eigenvalue of the unnormalized Hamiltonian

Notice: 
From the maximum eigenvalue of the Hamiltonian eN (in 'maxener') and the solution f(x)=e1 (in 'solution') of the optimization problem,
you can easily normalize the Hamiltonian, by dividing by (eN - e1)