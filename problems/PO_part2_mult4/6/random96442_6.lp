\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: - 698 x0 - 349 x1 - 258 x2 - 129 x3 - 96 x4 - 48 x5 + [ 472 x0^2
      + 472 x0*x1 + 384 x0*x2 + 192 x0*x3 + 176 x0*x4 + 88 x0*x5 + 118 x1^2
      + 192 x1*x2 + 96 x1*x3 + 88 x1*x4 + 44 x1*x5 + 808 x2^2 + 808 x2*x3
      + 96 x2*x4 + 48 x2*x5 + 202 x3^2 + 48 x3*x4 + 24 x3*x5 + 512 x4^2
      + 512 x4*x5 + 128 x5^2 ]/2
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
