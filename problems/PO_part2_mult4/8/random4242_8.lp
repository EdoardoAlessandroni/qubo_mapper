\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: - 228 x0 - 114 x1 - 312 x2 - 156 x3 - 288 x4 - 144 x5 - 262 x6 - 131 x7 +
      [ 632 x0^2 + 632 x0*x1 - 128 x0*x2 - 64 x0*x3 + 80 x0*x4 + 40 x0*x5
      + 1088 x0*x6 + 544 x0*x7 + 158 x1^2 - 64 x1*x2 - 32 x1*x3 + 40 x1*x4
      + 20 x1*x5 + 544 x1*x6 + 272 x1*x7 + 1168 x2^2 + 1168 x2*x3 + 80 x2*x4
      + 40 x2*x5 + 240 x2*x6 + 120 x2*x7 + 292 x3^2 + 40 x3*x4 + 20 x3*x5
      + 120 x3*x6 + 60 x3*x7 + 320 x4^2 + 320 x4*x5 + 128 x4*x6 + 64 x4*x7
      + 80 x5^2 + 64 x5*x6 + 32 x5*x7 + 696 x6^2 + 696 x6*x7 + 174 x7^2 ]/2
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
