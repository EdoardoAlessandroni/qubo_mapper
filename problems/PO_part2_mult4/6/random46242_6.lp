\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: - 340 x0 - 170 x1 - 420 x2 - 210 x3 + 356 x4 + 178 x5 + [ 528 x0^2
      + 528 x0*x1 + 64 x0*x2 + 32 x0*x3 + 384 x0*x4 + 192 x0*x5 + 132 x1^2
      + 32 x1*x2 + 16 x1*x3 + 192 x1*x4 + 96 x1*x5 + 376 x2^2 + 376 x2*x3
      + 544 x2*x4 + 272 x2*x5 + 94 x3^2 + 272 x3*x4 + 136 x3*x5 + 2656 x4^2
      + 2656 x4*x5 + 664 x5^2 ]/2
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
