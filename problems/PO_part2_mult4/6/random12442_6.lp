\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: 178 x0 + 89 x1 - 732 x2 - 366 x3 - 420 x4 - 210 x5 + [ 1064 x0^2
      + 1064 x0*x1 + 768 x0*x2 + 384 x0*x3 + 400 x0*x4 + 200 x0*x5 + 266 x1^2
      + 384 x1*x2 + 192 x1*x3 + 200 x1*x4 + 100 x1*x5 + 632 x2^2 + 632 x2*x3
      + 640 x2*x4 + 320 x2*x5 + 158 x3^2 + 320 x3*x4 + 160 x3*x5 + 376 x4^2
      + 376 x4*x5 + 94 x5^2 ]/2
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
