\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: - 258 x0 - 129 x1 + 628 x2 + 314 x3 - 852 x4 - 426 x5 + [ 808 x0^2
      + 808 x0*x1 - 48 x0*x2 - 24 x0*x3 + 208 x0*x4 + 104 x0*x5 + 202 x1^2
      - 24 x1*x2 - 12 x1*x3 + 104 x1*x4 + 52 x1*x5 + 2112 x2^2 + 2112 x2*x3
      + 1264 x2*x4 + 632 x2*x5 + 528 x3^2 + 632 x3*x4 + 316 x3*x5 + 1320 x4^2
      + 1320 x4*x5 + 330 x5^2 ]/2
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
