\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: - 156 x0 - 78 x1 - 290 x2 - 145 x3 + [ 504 x0^2 + 504 x0*x1 + 480 x0*x2
      + 240 x0*x3 + 126 x1^2 + 240 x1*x2 + 120 x1*x3 + 752 x2^2 + 752 x2*x3
      + 188 x3^2 ]/2
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
