\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: - 258 x0 - 129 x1 - 422 x2 - 211 x3 - 138 x4 - 69 x5 + [ 808 x0^2
      + 808 x0*x1 + 304 x0*x2 + 152 x0*x3 - 112 x0*x4 - 56 x0*x5 + 202 x1^2
      + 152 x1*x2 + 76 x1*x3 - 56 x1*x4 - 28 x1*x5 + 576 x2^2 + 576 x2*x3
      + 368 x2*x4 + 184 x2*x5 + 144 x3^2 + 184 x3*x4 + 92 x3*x5 + 592 x4^2
      + 592 x4*x5 + 148 x5^2 ]/2
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
