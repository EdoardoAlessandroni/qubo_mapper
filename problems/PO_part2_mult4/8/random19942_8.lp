\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: - 616 x0 - 308 x1 + 308 x2 + 154 x3 - 114 x4 - 57 x5 - 282 x6 - 141 x7 + [
      992 x0^2 + 992 x0*x1 + 880 x0*x2 + 440 x0*x3 + 736 x0*x6 + 368 x0*x7
      + 248 x1^2 + 440 x1*x2 + 220 x1*x3 + 368 x1*x6 + 184 x1*x7 + 1128 x2^2
      + 1128 x2*x3 + 144 x2*x4 + 72 x2*x5 + 864 x2*x6 + 432 x2*x7 + 282 x3^2
      + 72 x3*x4 + 36 x3*x5 + 432 x3*x6 + 216 x3*x7 + 328 x4^2 + 328 x4*x5
      + 304 x4*x6 + 152 x4*x7 + 82 x5^2 + 152 x5*x6 + 76 x5*x7 + 872 x6^2
      + 872 x6*x7 + 218 x7^2 ]/2
Subject To
 c0: 2 x0 + x1 + 2 x2 + x3 + 2 x4 + x5 + 2 x6 + x7 = 3

Bounds
 0 <= x0 <= 1
 0 <= x1 <= 1
 0 <= x2 <= 1
 0 <= x3 <= 1
 0 <= x4 <= 1
 0 <= x5 <= 1
 0 <= x6 <= 1
 0 <= x7 <= 1

Binaries
 x0 x1 x2 x3 x4 x5 x6 x7
End
