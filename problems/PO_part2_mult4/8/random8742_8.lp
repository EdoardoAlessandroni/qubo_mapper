\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: - 590 x0 - 295 x1 - 110 x2 - 55 x3 - 312 x4 - 156 x5 + 110 x6 + 55 x7 + [
      392 x0^2 + 392 x0*x1 + 1120 x0*x2 + 560 x0*x3 + 80 x0*x4 + 40 x0*x5
      + 432 x0*x6 + 216 x0*x7 + 98 x1^2 + 560 x1*x2 + 280 x1*x3 + 40 x1*x4
      + 20 x1*x5 + 216 x1*x6 + 108 x1*x7 + 1224 x2^2 + 1224 x2*x3 + 304 x2*x4
      + 152 x2*x5 + 688 x2*x6 + 344 x2*x7 + 306 x3^2 + 152 x3*x4 + 76 x3*x5
      + 344 x3*x6 + 172 x3*x7 + 1168 x4^2 + 1168 x4*x5 - 576 x4*x6 - 288 x4*x7
      + 292 x5^2 - 288 x5*x6 - 144 x5*x7 + 1192 x6^2 + 1192 x6*x7 + 298 x7^2 ]/2
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
