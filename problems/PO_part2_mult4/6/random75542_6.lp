\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: - 282 x0 - 141 x1 + 110 x2 + 55 x3 - 290 x4 - 145 x5 + [ 872 x0^2
      + 872 x0*x1 + 1632 x0*x2 + 816 x0*x3 + 1536 x0*x4 + 768 x0*x5 + 218 x1^2
      + 816 x1*x2 + 408 x1*x3 + 768 x1*x4 + 384 x1*x5 + 1192 x2^2 + 1192 x2*x3
      + 1520 x2*x4 + 760 x2*x5 + 298 x3^2 + 760 x3*x4 + 380 x3*x5 + 752 x4^2
      + 752 x4*x5 + 188 x5^2 ]/2
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
