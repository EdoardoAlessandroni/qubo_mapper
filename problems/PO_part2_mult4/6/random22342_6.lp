\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: - 434 x0 - 217 x1 + 104 x2 + 52 x3 + 586 x4 + 293 x5 + [ 3328 x0^2
      + 3328 x0*x1 + 1264 x0*x2 + 632 x0*x3 + 1456 x0*x4 + 728 x0*x5 + 832 x1^2
      + 632 x1*x2 + 316 x1*x3 + 728 x1*x4 + 364 x1*x5 + 2416 x2^2 + 2416 x2*x3
      - 320 x2*x4 - 160 x2*x5 + 604 x3^2 - 160 x3*x4 - 80 x3*x5 + 1240 x4^2
      + 1240 x4*x5 + 310 x5^2 ]/2
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
