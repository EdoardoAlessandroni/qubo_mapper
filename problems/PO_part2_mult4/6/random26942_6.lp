\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: 196 x0 + 98 x1 + 56 x2 + 28 x3 - 732 x4 - 366 x5 + [ 1216 x0^2
      + 1216 x0*x1 + 816 x0*x2 + 408 x0*x3 + 672 x0*x4 + 336 x0*x5 + 304 x1^2
      + 408 x1*x2 + 204 x1*x3 + 336 x1*x4 + 168 x1*x5 + 488 x2^2 + 488 x2*x3
      + 464 x2*x4 + 232 x2*x5 + 122 x3^2 + 232 x3*x4 + 116 x3*x5 + 632 x4^2
      + 632 x4*x5 + 158 x5^2 ]/2
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
