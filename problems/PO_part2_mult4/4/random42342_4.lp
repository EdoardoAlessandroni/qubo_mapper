\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: - 228 x0 - 114 x1 - 402 x2 - 201 x3 + [ 632 x0^2 + 632 x0*x1 + 560 x0*x2
      + 280 x0*x3 + 158 x1^2 + 280 x1*x2 + 140 x1*x3 + 424 x2^2 + 424 x2*x3
      + 106 x3^2 ]/2
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
