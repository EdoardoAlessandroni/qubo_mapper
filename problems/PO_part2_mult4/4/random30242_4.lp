\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: 172 x0 + 86 x1 - 422 x2 - 211 x3 + [ 1072 x0^2 + 1072 x0*x1 + 624 x0*x2
      + 312 x0*x3 + 268 x1^2 + 312 x1*x2 + 156 x1*x3 + 576 x2^2 + 576 x2*x3
      + 144 x3^2 ]/2
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
