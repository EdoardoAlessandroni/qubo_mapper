\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: - 1014 x0 - 507 x1 + 36 x2 + 18 x3 - 76 x4 - 38 x5 + [ 4880 x0^2
      + 4880 x0*x1 + 1376 x0*x2 + 688 x0*x3 + 480 x0*x4 + 240 x0*x5 + 1220 x1^2
      + 688 x1*x2 + 344 x1*x3 + 240 x1*x4 + 120 x1*x5 + 712 x2^2 + 712 x2*x3
      + 560 x2*x4 + 280 x2*x5 + 178 x3^2 + 280 x3*x4 + 140 x3*x5 + 560 x4^2
      + 560 x4*x5 + 140 x5^2 ]/2
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
