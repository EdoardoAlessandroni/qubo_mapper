\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: - 732 x0 - 366 x1 + 42 x2 + 21 x3 - 108 x4 - 54 x5 + [ 632 x0^2
      + 632 x0*x1 + 672 x0*x2 + 336 x0*x3 + 272 x0*x4 + 136 x0*x5 + 158 x1^2
      + 336 x1*x2 + 168 x1*x3 + 136 x1*x4 + 68 x1*x5 + 1144 x2^2 + 1144 x2*x3
      + 288 x2*x4 + 144 x2*x5 + 286 x3^2 + 144 x3*x4 + 72 x3*x5 + 360 x4^2
      + 360 x4*x5 + 90 x5^2 ]/2
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
