\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: - 672 x0 - 336 x1 + 374 x2 + 187 x3 + [ 1048 x0^2 + 1048 x0*x1 + 720 x0*x2
      + 360 x0*x3 + 262 x1^2 + 360 x1*x2 + 180 x1*x3 + 600 x2^2 + 600 x2*x3
      + 150 x3^2 ]/2
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
