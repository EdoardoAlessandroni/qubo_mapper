\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: - 68 x0 - 34 x1 - 66 x2 - 33 x3 + [ 2160 x0^2 + 2160 x0*x1 + 144 x0*x2
      + 72 x0*x3 + 540 x1^2 + 72 x1*x2 + 36 x1*x3 + 296 x2^2 + 296 x2*x3
      + 74 x3^2 ]/2
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
