\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: 18 x0 + 9 x1 + 562 x2 + 281 x3 + 196 x4 + 98 x5 + [ 832 x0^2 + 832 x0*x1
      + 368 x0*x2 + 184 x0*x3 + 1072 x0*x4 + 536 x0*x5 + 208 x1^2 + 184 x1*x2
      + 92 x1*x3 + 536 x1*x4 + 268 x1*x5 + 608 x2^2 + 608 x2*x3 + 208 x2*x4
      + 104 x2*x5 + 152 x3^2 + 104 x3*x4 + 52 x3*x5 + 1216 x4^2 + 1216 x4*x5
      + 304 x5^2 ]/2
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
