\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: - 114 x0 - 57 x1 + 208 x2 + 104 x3 + 200 x4 + 100 x5 + [ 328 x0^2
      + 328 x0*x1 + 32 x0*x2 + 16 x0*x3 + 256 x0*x4 + 128 x0*x5 + 82 x1^2
      + 16 x1*x2 + 8 x1*x3 + 128 x1*x4 + 64 x1*x5 + 680 x2^2 + 680 x2*x3
      + 352 x2*x4 + 176 x2*x5 + 170 x3^2 + 176 x3*x4 + 88 x3*x5 + 808 x4^2
      + 808 x4*x5 + 202 x5^2 ]/2
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
