\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: 40 x0 + 20 x1 + 200 x2 + 100 x3 - 68 x4 - 34 x5 + [ 984 x0^2 + 984 x0*x1
      + 592 x0*x2 + 296 x0*x3 + 1584 x0*x4 + 792 x0*x5 + 246 x1^2 + 296 x1*x2
      + 148 x1*x3 + 792 x1*x4 + 396 x1*x5 + 808 x2^2 + 808 x2*x3 + 528 x2*x4
      + 264 x2*x5 + 202 x3^2 + 264 x3*x4 + 132 x3*x5 + 2160 x4^2 + 2160 x4*x5
      + 540 x5^2 ]/2
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
