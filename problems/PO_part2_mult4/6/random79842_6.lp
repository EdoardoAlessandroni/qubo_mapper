\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: 178 x0 + 89 x1 + 18 x2 + 9 x3 - 168 x4 - 84 x5 + [ 1064 x0^2 + 1064 x0*x1
      + 656 x0*x2 + 328 x0*x3 + 368 x0*x4 + 184 x0*x5 + 266 x1^2 + 328 x1*x2
      + 164 x1*x3 + 184 x1*x4 + 92 x1*x5 + 832 x2^2 + 832 x2*x3 + 432 x2*x4
      + 216 x2*x5 + 208 x3^2 + 216 x3*x4 + 108 x3*x5 + 456 x4^2 + 456 x4*x5
      + 114 x5^2 ]/2
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
