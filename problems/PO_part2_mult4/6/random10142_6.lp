\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: - 672 x0 - 336 x1 + 52 x2 + 26 x3 - 254 x4 - 127 x5 + [ 1048 x0^2
      + 1048 x0*x1 + 768 x0*x2 + 384 x0*x3 + 272 x0*x4 + 136 x0*x5 + 262 x1^2
      + 384 x1*x2 + 192 x1*x3 + 136 x1*x4 + 68 x1*x5 + 2448 x2^2 + 2448 x2*x3
      + 800 x2*x4 + 400 x2*x5 + 612 x3^2 + 400 x3*x4 + 200 x3*x5 + 176 x4^2
      + 176 x4*x5 + 44 x5^2 ]/2
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
