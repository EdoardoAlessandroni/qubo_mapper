\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: 374 x0 + 187 x1 + 196 x2 + 98 x3 + 214 x4 + 107 x5 + 56 x6 + 28 x7 + [
      600 x0^2 + 600 x0*x1 + 992 x0*x2 + 496 x0*x3 + 832 x0*x4 + 416 x0*x5
      + 448 x0*x6 + 224 x0*x7 + 150 x1^2 + 496 x1*x2 + 248 x1*x3 + 416 x1*x4
      + 208 x1*x5 + 224 x1*x6 + 112 x1*x7 + 1216 x2^2 + 1216 x2*x3 + 1024 x2*x4
      + 512 x2*x5 + 816 x2*x6 + 408 x2*x7 + 304 x3^2 + 512 x3*x4 + 256 x3*x5
      + 408 x3*x6 + 204 x3*x7 + 1104 x4^2 + 1104 x4*x5 + 608 x4*x6 + 304 x4*x7
      + 276 x5^2 + 304 x5*x6 + 152 x5*x7 + 488 x6^2 + 488 x6*x7 + 122 x7^2 ]/2
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
