\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: 66 x0 + 33 x1 - 152 x2 - 76 x3 + [ 944 x0^2 + 944 x0*x1 + 336 x0*x2
      + 168 x0*x3 + 236 x1^2 + 168 x1*x2 + 84 x1*x3 + 480 x2^2 + 480 x2*x3
      + 120 x3^2 ]/2
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
