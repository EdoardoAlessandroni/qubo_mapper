\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: - 590 x0 - 295 x1 + 40 x2 + 20 x3 - 104 x4 - 52 x5 - 236 x6 - 118 x7 + [
      392 x0^2 + 392 x0*x1 + 768 x0*x2 + 384 x0*x3 + 720 x0*x4 + 360 x0*x5
      + 352 x0*x6 + 176 x0*x7 + 98 x1^2 + 384 x1*x2 + 192 x1*x3 + 360 x1*x4
      + 180 x1*x5 + 176 x1*x6 + 88 x1*x7 + 984 x2^2 + 984 x2*x3 + 1424 x2*x4
      + 712 x2*x5 + 1840 x2*x6 + 920 x2*x7 + 246 x3^2 + 712 x3*x4 + 356 x3*x5
      + 920 x3*x6 + 460 x3*x7 + 2608 x4^2 + 2608 x4*x5 + 560 x4*x6 + 280 x4*x7
      + 652 x5^2 + 280 x5*x6 + 140 x5*x7 + 2448 x6^2 + 2448 x6*x7 + 612 x7^2 ]/2
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
