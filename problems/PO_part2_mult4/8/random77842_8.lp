\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: 400 x0 + 200 x1 - 434 x2 - 217 x3 + 72 x4 + 36 x5 + 42 x6 + 21 x7 + [
      2528 x0^2 + 2528 x0*x1 + 768 x0*x2 + 384 x0*x3 + 48 x0*x4 + 24 x0*x5
      + 1520 x0*x6 + 760 x0*x7 + 632 x1^2 + 384 x1*x2 + 192 x1*x3 + 24 x1*x4
      + 12 x1*x5 + 760 x1*x6 + 380 x1*x7 + 3328 x2^2 + 3328 x2*x3 + 896 x2*x4
      + 448 x2*x5 + 1120 x2*x6 + 560 x2*x7 + 832 x3^2 + 448 x3*x4 + 224 x3*x5
      + 560 x3*x6 + 280 x3*x7 + 552 x4^2 + 552 x4*x5 + 736 x4*x6 + 368 x4*x7
      + 138 x5^2 + 368 x5*x6 + 184 x5*x7 + 1016 x6^2 + 1016 x6*x7 + 254 x7^2 ]/2
Subject To
 c0: 2 x0 + x1 + 2 x2 + x3 + 2 x4 + x5 + 2 x6 + x7 = 3

Bounds
 0 <= x0 <= 1
 0 <= x1 <= 1
 0 <= x2 <= 1
 0 <= x3 <= 1
 0 <= x4 <= 1
 0 <= x5 <= 1
 0 <= x6 <= 1
 0 <= x7 <= 1

Binaries
 x0 x1 x2 x3 x4 x5 x6 x7
End
