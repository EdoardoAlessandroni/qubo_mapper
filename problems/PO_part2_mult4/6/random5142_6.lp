\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: 90 x0 + 45 x1 + 90 x2 + 45 x3 + 196 x4 + 98 x5 + [ 664 x0^2 + 664 x0*x1
      + 784 x0*x2 + 392 x0*x3 + 544 x0*x4 + 272 x0*x5 + 166 x1^2 + 392 x1*x2
      + 196 x1*x3 + 272 x1*x4 + 136 x1*x5 + 448 x2^2 + 448 x2*x3 + 672 x2*x4
      + 336 x2*x5 + 112 x3^2 + 336 x3*x4 + 168 x3*x5 + 1216 x4^2 + 1216 x4*x5
      + 304 x5^2 ]/2
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
