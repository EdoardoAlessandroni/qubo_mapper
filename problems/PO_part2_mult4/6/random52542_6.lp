\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: - 418 x0 - 209 x1 - 66 x2 - 33 x3 - 90 x4 - 45 x5 + [ 2128 x0^2
      + 2128 x0*x1 - 288 x0*x2 - 144 x0*x3 + 1152 x0*x4 + 576 x0*x5 + 532 x1^2
      - 144 x1*x2 - 72 x1*x3 + 576 x1*x4 + 288 x1*x5 + 296 x2^2 + 296 x2*x3
      + 96 x2*x4 + 48 x2*x5 + 74 x3^2 + 48 x3*x4 + 24 x3*x5 + 664 x4^2
      + 664 x4*x5 + 166 x5^2 ]/2
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
