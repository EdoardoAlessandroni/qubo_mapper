\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random1742_9_1.lp

Minimize
 obj: 6 x0 + 4 x1 + 3 x2 + 8 x3 + 10 x4 - 2 x5 - 2 x6 - 4 x7 - 7 x8 + [ 2 x0*x1
      - 8 x0*x2 - 10 x0*x4 + 10 x0*x6 - 18 x0*x8 + 4 x1*x2 + 20 x1*x3 + 2 x1*x4
      - 18 x1*x5 - 12 x2*x3 + 10 x2*x5 + 12 x2*x7 - 14 x3*x4 + 8 x3*x5
      + 16 x3*x7 - 12 x4*x6 - 8 x4*x8 + 6 x5*x7 - 20 x5*x8 - 10 x6*x7 - 10 x6*x8
      - 16 x7*x8 ]/2
Subject To
 c0: x1 - x2 - x3 + 2 x4 - x8 = 0

Bounds
 0 <= x0 <= 1
 0 <= x1 <= 1
 0 <= x2 <= 1
 0 <= x3 <= 1
 0 <= x4 <= 1
 0 <= x5 <= 1
 0 <= x6 <= 1
 0 <= x7 <= 1
 0 <= x8 <= 1

Binaries
 x0 x1 x2 x3 x4 x5 x6 x7 x8
End
