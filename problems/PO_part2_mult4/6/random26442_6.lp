\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: - 286 x0 - 143 x1 + 104 x2 + 52 x3 - 200 x4 - 100 x5 + [ 728 x0^2
      + 728 x0*x1 - 16 x0*x2 - 8 x0*x3 - 48 x0*x4 - 24 x0*x5 + 182 x1^2
      - 8 x1*x2 - 4 x1*x3 - 24 x1*x4 - 12 x1*x5 + 2416 x2^2 + 2416 x2*x3
      + 912 x2*x4 + 456 x2*x5 + 604 x3^2 + 456 x3*x4 + 228 x3*x5 + 248 x4^2
      + 248 x4*x5 + 62 x5^2 ]/2
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
