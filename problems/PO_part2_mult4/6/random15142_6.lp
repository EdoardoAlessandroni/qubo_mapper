\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: 252 x0 + 126 x1 + 284 x2 + 142 x3 - 156 x4 - 78 x5 + [ 536 x0^2
      + 536 x0*x1 + 416 x0*x2 + 208 x0*x3 + 400 x0*x4 + 200 x0*x5 + 134 x1^2
      + 208 x1*x2 + 104 x1*x3 + 200 x1*x4 + 100 x1*x5 + 1104 x2^2 + 1104 x2*x3
      + 768 x2*x4 + 384 x2*x5 + 276 x3^2 + 384 x3*x4 + 192 x3*x5 + 504 x4^2
      + 504 x4*x5 + 126 x5^2 ]/2
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
