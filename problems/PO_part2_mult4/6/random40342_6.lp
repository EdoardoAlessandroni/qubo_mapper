\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: - 110 x0 - 55 x1 - 196 x2 - 98 x3 + 370 x4 + 185 x5 + [ 1224 x0^2
      + 1224 x0*x1 + 944 x0*x2 + 472 x0*x3 + 1344 x0*x4 + 672 x0*x5 + 306 x1^2
      + 472 x1*x2 + 236 x1*x3 + 672 x1*x4 + 336 x1*x5 + 344 x2^2 + 344 x2*x3
      + 496 x2*x4 + 248 x2*x5 + 86 x3^2 + 248 x3*x4 + 124 x3*x5 + 960 x4^2
      + 960 x4*x5 + 240 x5^2 ]/2
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
