\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: 42 x0 + 21 x1 + 586 x2 + 293 x3 - 52 x4 - 26 x5 + [ 1144 x0^2 + 1144 x0*x1
      - 480 x0*x2 - 240 x0*x3 + 224 x0*x4 + 112 x0*x5 + 286 x1^2 - 240 x1*x2
      - 120 x1*x3 + 112 x1*x4 + 56 x1*x5 + 1240 x2^2 + 1240 x2*x3 + 336 x2*x4
      + 168 x2*x5 + 310 x3^2 + 168 x3*x4 + 84 x3*x5 + 376 x4^2 + 376 x4*x5
      + 94 x5^2 ]/2
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
