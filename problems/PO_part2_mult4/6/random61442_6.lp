\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: - 590 x0 - 295 x1 - 298 x2 - 149 x3 - 296 x4 - 148 x5 + [ 392 x0^2
      + 392 x0*x1 + 528 x0*x2 + 264 x0*x3 + 544 x0*x4 + 272 x0*x5 + 98 x1^2
      + 264 x1*x2 + 132 x1*x3 + 272 x1*x4 + 136 x1*x5 + 768 x2^2 + 768 x2*x3
      + 352 x2*x4 + 176 x2*x5 + 192 x3^2 + 176 x3*x4 + 88 x3*x5 + 536 x4^2
      + 536 x4*x5 + 134 x5^2 ]/2
Subject To
 c0: 2 x0 + x1 + 2 x2 + x3 + 2 x4 + x5 = 3

Bounds
 0 <= x0 <= 1
 0 <= x1 <= 1
 0 <= x2 <= 1
 0 <= x3 <= 1
 0 <= x4 <= 1
 0 <= x5 <= 1

Binaries
 x0 x1 x2 x3 x4 x5
End
