\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: 376 x0 + 188 x1 - 856 x2 - 428 x3 - 288 x4 - 144 x5 + [ 792 x0^2
      + 792 x0*x1 + 1792 x0*x2 + 896 x0*x3 - 176 x0*x4 - 88 x0*x5 + 198 x1^2
      + 896 x1*x2 + 448 x1*x3 - 88 x1*x4 - 44 x1*x5 + 3232 x2^2 + 3232 x2*x3
      - 432 x2*x4 - 216 x2*x5 + 808 x3^2 - 216 x3*x4 - 108 x3*x5 + 320 x4^2
      + 320 x4*x5 + 80 x5^2 ]/2
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
