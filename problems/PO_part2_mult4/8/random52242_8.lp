\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: 10 x0 + 5 x1 - 88 x2 - 44 x3 - 200 x4 - 100 x5 + 196 x6 + 98 x7 + [
      1000 x0^2 + 1000 x0*x1 + 480 x0*x2 + 240 x0*x3 + 448 x0*x4 + 224 x0*x5
      + 1408 x0*x6 + 704 x0*x7 + 250 x1^2 + 240 x1*x2 + 120 x1*x3 + 224 x1*x4
      + 112 x1*x5 + 704 x1*x6 + 352 x1*x7 + 944 x2^2 + 944 x2*x3 + 224 x2*x4
      + 112 x2*x5 + 752 x2*x6 + 376 x2*x7 + 236 x3^2 + 112 x3*x4 + 56 x3*x5
      + 376 x3*x6 + 188 x3*x7 + 248 x4^2 + 248 x4*x5 + 176 x4*x6 + 88 x4*x7
      + 62 x5^2 + 88 x5*x6 + 44 x5*x7 + 1216 x6^2 + 1216 x6*x7 + 304 x7^2 ]/2
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
