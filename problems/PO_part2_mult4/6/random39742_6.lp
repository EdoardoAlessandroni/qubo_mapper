\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: - 434 x0 - 217 x1 - 88 x2 - 44 x3 + 10 x4 + 5 x5 + [ 3328 x0^2
      + 3328 x0*x1 + 768 x0*x2 + 384 x0*x3 + 256 x0*x4 + 128 x0*x5 + 832 x1^2
      + 384 x1*x2 + 192 x1*x3 + 128 x1*x4 + 64 x1*x5 + 472 x2^2 + 472 x2*x3
      + 560 x2*x4 + 280 x2*x5 + 118 x3^2 + 280 x3*x4 + 140 x3*x5 + 1000 x4^2
      + 1000 x4*x5 + 250 x5^2 ]/2
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
