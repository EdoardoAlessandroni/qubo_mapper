\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: - 78 x0 - 39 x1 + 66 x2 + 33 x3 + 42 x4 + 21 x5 + [ 432 x0^2 + 432 x0*x1
      + 656 x0*x2 + 328 x0*x3 + 912 x0*x4 + 456 x0*x5 + 108 x1^2 + 328 x1*x2
      + 164 x1*x3 + 456 x1*x4 + 228 x1*x5 + 944 x2^2 + 944 x2*x3 + 960 x2*x4
      + 480 x2*x5 + 236 x3^2 + 480 x3*x4 + 240 x3*x5 + 1144 x4^2 + 1144 x4*x5
      + 286 x5^2 ]/2
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
