\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: 216 x0 + 108 x1 - 458 x2 - 229 x3 + 214 x4 + 107 x5 + [ 1144 x0^2
      + 1144 x0*x1 + 400 x0*x2 + 200 x0*x3 + 992 x0*x4 + 496 x0*x5 + 286 x1^2
      + 200 x1*x2 + 100 x1*x3 + 496 x1*x4 + 248 x1*x5 + 816 x2^2 + 816 x2*x3
      - 160 x2*x4 - 80 x2*x5 + 204 x3^2 - 80 x3*x4 - 40 x3*x5 + 1104 x4^2
      + 1104 x4*x5 + 276 x5^2 ]/2
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
