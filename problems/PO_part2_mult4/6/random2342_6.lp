\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: - 262 x0 - 131 x1 - 290 x2 - 145 x3 - 1014 x4 - 507 x5 + [ 696 x0^2
      + 696 x0*x1 + 672 x0*x2 + 336 x0*x3 + 1792 x0*x4 + 896 x0*x5 + 174 x1^2
      + 336 x1*x2 + 168 x1*x3 + 896 x1*x4 + 448 x1*x5 + 752 x2^2 + 752 x2*x3
      - 80 x2*x4 - 40 x2*x5 + 188 x3^2 - 40 x3*x4 - 20 x3*x5 + 4880 x4^2
      + 4880 x4*x5 + 1220 x5^2 ]/2
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
