Samples in folders are named as random"seed"_"nvars"_"nconstraints".lp, for example random42_50_12.lp



FOLDERS


################ NN ################


NN_linear_deg5
Random sparse LCBO (NN stands for No Normalization: objective function stays with interger coefficients)
Bvars:	range(4, 26)
N_constraints:	max(1, int(n_vars/5))
Densities:	inversely proportional to n (5/n in particular), since every node has a maximum degree and thus the number of terms in the objective function is linear (from which the folder name) with n.
Maximum degree: 5
Number of samples = 1000

NN_linear_deg10
Random sparse LCBO (NN stands for No Normalization: objective function stays with interger coefficients)
Bvars:	range(4, 26) + range(30, 76, 5)
N_constraints:	max(1, int(n_vars/5))
Densities:	inversely proportional to n (10/n in particular), since every node has a maximum degree and thus the number of terms in the objective function is linear (from which the folder name) with n.
Maximum degree: 10
Number of samples = 1000


################ SPP ################


SPP_p15.txt
Set Partitioning Problem instances
Bvars: range(6,31)
N_constraints:	max(1, int(n_vars/3))
1000 samples
Density: obj is linear (vector full). Constraint matrix A is sampled with density (prob A_ij = 1) = 0.15

SPP_p25.txt
Set Partitioning Problem instances
Bvars: range(6,31)
N_constraints:	max(1, int(n_vars/3))
1000 samples
Density: obj is linear (vector full). Constraint matrix A is sampled with density (prob A_ij = 1) = 0.25

SPP_p35.txt
Set Partitioning Problem instances
Bvars: range(6,31)
N_constraints:	max(1, int(n_vars/3))
1000 samples
Density: obj is linear (vector full). Constraint matrix A is sampled with density (prob A_ij = 1) = 0.35
SPP stands for Set Partitioning Problem


################ PO ################


PO_part*
Portfolio Optimization, Markowitz model. Data from S&P500 financial data.
n = range(6,25,3) for w=3.	range(6,25,3) for w=2.
1000 samples
pars: w=[part](2, 3), multiplier=[mult](1e2, 1e4), risk aversion=[ra](0.5, 1, 2. If not specified, 1.)

PO_norm_part*
Portfolio Optimization, Markowitz model. Data from S&P500 financial data. , (NORMalized formulation, i.e. decision variables properly describe fractions)
n = range(6,25,3)
1000 samples
pars: w=[part](3), multiplier=[mult](1e2, 1e4), risk aversion=1 (since not specified)

PO_norm_big
Portfolio Optimization, Markowitz model. Data from S&P500 financial data. (NORMalized formulation). Big number of variables.
range(10,51,10) + [75] + range(100, 401, 50)
100 samples
pars: w=5, multiplier=1e4, risk aversion=1





################ Old (deleted) ################


NN_25 (deleted)
Bvars:	range(4,25,2) + [30, 35, 40, 45, 50]
N_constraints:	max(1, int(n_vars/4))
Densities:	.25 for Q, L in objective function and for A constraint matrix
Number of samples = 100
NN stands for No Normalization (objective function stays with interger coefficients)

NN_50 (deleted)
Bvars:	range(4,25,2) + [30, 35, 40, 45, 50]
N_constraints:	max(1, int(n_vars/4))
Densities:	.5 for Q, L in objective function and for A constraint matrix
Number of samples = 100
NN stands for No Normalization (objective function stays with interger coefficients)

NN (deleted)
Bvars:	range(2,25,2) + [30, 35, 40, 45, 50]
N_constraints:	max(1, int(n_vars/4))
Densities:	1 for Q, L in objective function, .5 for A constraint matrix
NN stands for No Normalization (objective function stays with interger coefficients)

SN (deleted)
Bvars:	range(2,21,2) + [24, 30, 35, 40, 45, 50]
N_constraints:	max(1, int(n_vars/4))
Densities:	1 for Q, L in objective function, .5 for A constraint matrix
Number of samples = 100
SN stands for Sampling Normalization (bounding spectral norm of Hamiltonian by sampling coefficients in objective function ~ 1/n^2)

BN (deleted)
Bvars:	range(2,29,2)
N_constraints:	max(1, int(n_vars/4)).
Densities:	1 for Q, L in objective function, .5 for A constraint matrix
Number of samples = 100
BN stands for Brute Normalization (computing full hamiltoinian and dividing by spectral gap)

old (deleted)
Bvars:	range(2,21)
N_constraints:	max(1, int(n_vars/3)).
Densities:	1 for Q, L in objective function, .5 for A constraint matrix
Number of samples = 100
Brute Normalization was emploied (computing full hamiltoinian and dividing by spectral gap)
