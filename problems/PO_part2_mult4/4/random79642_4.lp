\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: 562 x0 + 281 x1 - 262 x2 - 131 x3 + [ 608 x0^2 + 608 x0*x1 + 624 x0*x2
      + 312 x0*x3 + 152 x1^2 + 312 x1*x2 + 156 x1*x3 + 696 x2^2 + 696 x2*x3
      + 174 x3^2 ]/2
Subject To
 c0: 2 x0 + x1 + 2 x2 + x3 = 3

Bounds
 0 <= x0 <= 1
 0 <= x1 <= 1
 0 <= x2 <= 1
 0 <= x3 <= 1

Binaries
 x0 x1 x2 x3
End
