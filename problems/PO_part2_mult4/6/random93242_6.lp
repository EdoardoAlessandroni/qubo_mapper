\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: - 114 x0 - 57 x1 + 584 x2 + 292 x3 + 90 x4 + 45 x5 + [ 328 x0^2
      + 328 x0*x1 + 480 x0*x2 + 240 x0*x3 + 368 x0*x4 + 184 x0*x5 + 82 x1^2
      + 240 x1*x2 + 120 x1*x3 + 184 x1*x4 + 92 x1*x5 + 1304 x2^2 + 1304 x2*x3
      + 800 x2*x4 + 400 x2*x5 + 326 x3^2 + 400 x3*x4 + 200 x3*x5 + 2528 x4^2
      + 2528 x4*x5 + 632 x5^2 ]/2
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
