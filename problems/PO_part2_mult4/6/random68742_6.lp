\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: 328 x0 + 164 x1 - 156 x2 - 78 x3 - 66 x4 - 33 x5 + [ 1128 x0^2
      + 1128 x0*x1 + 816 x0*x2 + 408 x0*x3 + 640 x0*x4 + 320 x0*x5 + 282 x1^2
      + 408 x1*x2 + 204 x1*x3 + 320 x1*x4 + 160 x1*x5 + 504 x2^2 + 504 x2*x3
      + 272 x2*x4 + 136 x2*x5 + 126 x3^2 + 136 x3*x4 + 68 x3*x5 + 296 x4^2
      + 296 x4*x5 + 74 x5^2 ]/2
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
