\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: - 78 x0 - 39 x1 - 672 x2 - 336 x3 + 870 x4 + 435 x5 - 616 x6 - 308 x7 + [
      432 x0^2 + 432 x0*x1 + 560 x0*x2 + 280 x0*x3 + 784 x0*x4 + 392 x0*x5
      + 544 x0*x6 + 272 x0*x7 + 108 x1^2 + 280 x1*x2 + 140 x1*x3 + 392 x1*x4
      + 196 x1*x5 + 272 x1*x6 + 136 x1*x7 + 1048 x2^2 + 1048 x2*x3 + 512 x2*x4
      + 256 x2*x5 + 624 x2*x6 + 312 x2*x7 + 262 x3^2 + 256 x3*x4 + 128 x3*x5
      + 312 x3*x6 + 156 x3*x7 + 1032 x4^2 + 1032 x4*x5 + 832 x4*x6 + 416 x4*x7
      + 258 x5^2 + 416 x5*x6 + 208 x5*x7 + 992 x6^2 + 992 x6*x7 + 248 x7^2 ]/2
Subject To
 c0: 2 x0 + x1 + 2 x2 + x3 + 2 x4 + x5 + 2 x6 + x7 = 3

Bounds
 0 <= x0 <= 1
 0 <= x1 <= 1
 0 <= x2 <= 1
 0 <= x3 <= 1
 0 <= x4 <= 1
 0 <= x5 <= 1
 0 <= x6 <= 1
 0 <= x7 <= 1

Binaries
 x0 x1 x2 x3 x4 x5 x6 x7
End
