\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: 90 x0 + 45 x1 - 66 x2 - 33 x3 - 258 x4 - 129 x5 - 10 x6 - 5 x7 + [
      664 x0^2 + 664 x0*x1 + 176 x0*x2 + 88 x0*x3 - 208 x0*x4 - 104 x0*x5
      + 208 x0*x6 + 104 x0*x7 + 166 x1^2 + 88 x1*x2 + 44 x1*x3 - 104 x1*x4
      - 52 x1*x5 + 104 x1*x6 + 52 x1*x7 + 296 x2^2 + 296 x2*x3 + 144 x2*x4
      + 72 x2*x5 + 576 x2*x6 + 288 x2*x7 + 74 x3^2 + 72 x3*x4 + 36 x3*x5
      + 288 x3*x6 + 144 x3*x7 + 808 x4^2 + 808 x4*x5 + 256 x4*x6 + 128 x4*x7
      + 202 x5^2 + 128 x5*x6 + 64 x5*x7 + 912 x6^2 + 912 x6*x7 + 228 x7^2 ]/2
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
