\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: - 76 x0 - 38 x1 + 206 x2 + 103 x3 - 590 x4 - 295 x5 - 66 x6 - 33 x7 + [
      560 x0^2 + 560 x0*x1 + 496 x0*x2 + 248 x0*x3 + 400 x0*x4 + 200 x0*x5
      + 496 x0*x6 + 248 x0*x7 + 140 x1^2 + 248 x1*x2 + 124 x1*x3 + 200 x1*x4
      + 100 x1*x5 + 248 x1*x6 + 124 x1*x7 + 984 x2^2 + 984 x2*x3 + 848 x2*x4
      + 424 x2*x5 + 192 x2*x6 + 96 x2*x7 + 246 x3^2 + 424 x3*x4 + 212 x3*x5
      + 96 x3*x6 + 48 x3*x7 + 392 x4^2 + 392 x4*x5 + 48 x4*x6 + 24 x4*x7
      + 98 x5^2 + 24 x5*x6 + 12 x5*x7 + 296 x6^2 + 296 x6*x7 + 74 x7^2 ]/2
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
