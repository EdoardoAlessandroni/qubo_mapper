\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: - 852 x0 - 426 x1 - 88 x2 - 44 x3 + [ 1320 x0^2 + 1320 x0*x1 + 208 x0*x2
      + 104 x0*x3 + 330 x1^2 + 104 x1*x2 + 52 x1*x3 + 472 x2^2 + 472 x2*x3
      + 118 x3^2 ]/2
Subject To
 c0: 2 x0 + x1 + 2 x2 + x3 = 3

Bounds
 0 <= x0 <= 1
 0 <= x1 <= 1
 0 <= x2 <= 1
 0 <= x3 <= 1

Binaries
 x0 x1 x2 x3
End
