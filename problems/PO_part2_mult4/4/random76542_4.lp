\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: 72 x0 + 36 x1 + 200 x2 + 100 x3 + [ 552 x0^2 + 552 x0*x1 + 208 x0*x2
      + 104 x0*x3 + 138 x1^2 + 104 x1*x2 + 52 x1*x3 + 808 x2^2 + 808 x2*x3
      + 202 x3^2 ]/2
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
