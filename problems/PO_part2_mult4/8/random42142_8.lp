\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: 90 x0 + 45 x1 - 430 x2 - 215 x3 - 10 x4 - 5 x5 + 400 x6 + 200 x7 + [
      448 x0^2 + 448 x0*x1 + 208 x0*x2 + 104 x0*x3 + 608 x0*x4 + 304 x0*x5
      + 1120 x0*x6 + 560 x0*x7 + 112 x1^2 + 104 x1*x2 + 52 x1*x3 + 304 x1*x4
      + 152 x1*x5 + 560 x1*x6 + 280 x1*x7 + 512 x2^2 + 512 x2*x3 - 336 x2*x4
      - 168 x2*x5 - 160 x2*x6 - 80 x2*x7 + 128 x3^2 - 168 x3*x4 - 84 x3*x5
      - 80 x3*x6 - 40 x3*x7 + 912 x4^2 + 912 x4*x5 + 2736 x4*x6 + 1368 x4*x7
      + 228 x5^2 + 1368 x5*x6 + 684 x5*x7 + 2528 x6^2 + 2528 x6*x7 + 632 x7^2
      ]/2
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
