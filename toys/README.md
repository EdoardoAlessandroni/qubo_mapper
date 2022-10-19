Samples in folders are named as random"seed"_"nvars"_"nconstraints".lp, for example random42_50_12.lp



FOLDERS

SPP_deg10
Bvars:	range(6, 26)
N_constraints:	max(1, int(n_vars/3))
Densities:	Objective function is only linear, constraint matrix A is sampled using a probability p (\sim prob A_ij = 1) which is set as follows:
For n = 6,...,19 p = 0.35, and for this value the quadratic part in the penalized function A^tA is such that the resulting degree of the nodes is on average smaller than degree_max = 10.
For n = 20,...,25 p is given by the following formula, which ensures that the quadratic part in the penalized function A^tA is such that the resulting degree of the nodes is on average equal to degree_max = 10.
p(delta) = sqrt( 1 - power( 1 - delta/(n-1) , 1/m ) )
Where n = number of binary variables and m is the number of constraints m = int( n/3 )
For a better discussion about where this formula arises from, see the file candidate_problems.md.
Maximum average degree: 10
Number of samples = 1000
SPP stands for Set Partitioning Problem

SPP_deg5
Bvars:	range(6, 26)
N_constraints:	max(1, int(n_vars/3))
Densities:	Objective function is only linear, constraint matrix A is sampled using a probability p (\sim prob A_ij = 1) which is set as follows:
For n = 6,...,13 p = 0.35, and for this value the quadratic part in the penalized function A^tA is such that the resulting degree of the nodes is on average smaller than degree_max = 5.
For n = 14,...,25 p is given by the following formula, which ensures that the quadratic part in the penalized function A^tA is such that the resulting degree of the nodes is on average equal to degree_max = 5.
p(delta) = sqrt( 1 - power( 1 - delta/(n-1) , 1/m ) )
Where n = number of binary variables and m is the number of constraints m = int( n/3 )
For a better discussion about where this formula arises from, see the file candidate_problems.md.
Maximum average degree: 5
Number of samples = 1000
SPP stands for Set Partitioning Problem

NN_linear_deg5
Bvars:	range(4, 26)
N_constraints:	max(1, int(n_vars/5))
Densities:	not constant, but rather inversely proportional to n (5/n in particular). This is because every node has a maximum degree and thus the number of terms in the objective function is linear (from which the folder name) with n.
Maximum degree: 5
Number of samples = 1000
NN stands for No Normalization (objective function stays with interger coefficients)


NN_linear_deg10
Bvars:	range(4, 26) + range(30, 76, 5)
N_constraints:	max(1, int(n_vars/5))
Densities:	not constant, but rather inversely proportional to n (10/n in particular). This is because every node has a maximum degree and thus the number of terms in the objective function is linear (from which the folder name) with n.
Maximum degree: 10
Number of samples = 1000
NN stands for No Normalization (objective function stays with interger coefficients)


NN_25
Bvars:	range(4,25,2) + [30, 35, 40, 45, 50]
N_constraints:	max(1, int(n_vars/4))
Densities:	.25 for Q, L in objective function and for A constraint matrix
Number of samples = 100
NN stands for No Normalization (objective function stays with interger coefficients)


NN_50
Bvars:	range(4,25,2) + [30, 35, 40, 45, 50]
N_constraints:	max(1, int(n_vars/4))
Densities:	.5 for Q, L in objective function and for A constraint matrix
Number of samples = 100
NN stands for No Normalization (objective function stays with interger coefficients)


NN
Bvars:	range(2,25,2) + [30, 35, 40, 45, 50]
N_constraints:	max(1, int(n_vars/4))
Densities:	1 for Q, L in objective function, .5 for A constraint matrix
NN stands for No Normalization (objective function stays with interger coefficients)


SN
Bvars:	range(2,21,2) + [24, 30, 35, 40, 45, 50]
N_constraints:	max(1, int(n_vars/4))
Densities:	1 for Q, L in objective function, .5 for A constraint matrix
Number of samples = 100
SN stands for Sampling Normalization (bounding spectral norm of Hamiltonian by sampling coefficients in objective function ~ 1/n^2)


BN
Bvars:	range(2,29,2)
N_constraints:	max(1, int(n_vars/4)).
Densities:	1 for Q, L in objective function, .5 for A constraint matrix
Number of samples = 100
BN stands for Brute Normalization (computing full hamiltoinian and dividing by spectral gap)

old
Bvars:	range(2,21)
N_constraints:	max(1, int(n_vars/3)).
Densities:	1 for Q, L in objective function, .5 for A constraint matrix
Number of samples = 100
Brute Normalization was emploied (computing full hamiltoinian and dividing by spectral gap)
