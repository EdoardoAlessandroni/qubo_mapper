\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: 230 x0 + 115 x1 - 52 x2 - 26 x3 + 70 x4 + 35 x5 + [ 1000 x0^2 + 1000 x0*x1
      + 464 x0*x2 + 232 x0*x3 + 640 x0*x4 + 320 x0*x5 + 250 x1^2 + 232 x1*x2
      + 116 x1*x3 + 320 x1*x4 + 160 x1*x5 + 376 x2^2 + 376 x2*x3 + 400 x2*x4
      + 200 x2*x5 + 94 x3^2 + 200 x3*x4 + 100 x3*x5 + 440 x4^2 + 440 x4*x5
      + 110 x5^2 ]/2
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
