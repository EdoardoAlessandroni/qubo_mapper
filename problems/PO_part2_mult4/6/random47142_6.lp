\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: 278 x0 + 139 x1 + 370 x2 + 185 x3 - 236 x4 - 118 x5 + [ 600 x0^2
      + 600 x0*x1 + 512 x0*x2 + 256 x0*x3 + 160 x0*x4 + 80 x0*x5 + 150 x1^2
      + 256 x1*x2 + 128 x1*x3 + 80 x1*x4 + 40 x1*x5 + 960 x2^2 + 960 x2*x3
      + 224 x2*x4 + 112 x2*x5 + 240 x3^2 + 112 x3*x4 + 56 x3*x5 + 464 x4^2
      + 464 x4*x5 + 116 x5^2 ]/2
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
