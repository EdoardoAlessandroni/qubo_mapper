\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: - 200 x0 - 100 x1 - 672 x2 - 336 x3 - 78 x4 - 39 x5 + [ 248 x0^2
      + 248 x0*x1 + 32 x0*x2 + 16 x0*x3 + 256 x0*x4 + 128 x0*x5 + 62 x1^2
      + 16 x1*x2 + 8 x1*x3 + 128 x1*x4 + 64 x1*x5 + 1048 x2^2 + 1048 x2*x3
      + 560 x2*x4 + 280 x2*x5 + 262 x3^2 + 280 x3*x4 + 140 x3*x5 + 432 x4^2
      + 432 x4*x5 + 108 x5^2 ]/2
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
