\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: - 10 x0 - 5 x1 + 90 x2 + 45 x3 - 1014 x4 - 507 x5 + 356 x6 + 178 x7
      + 498 x8 + 249 x9 + [ 912 x0^2 + 912 x0*x1 + 608 x0*x2 + 304 x0*x3
      - 592 x0*x4 - 296 x0*x5 + 480 x0*x6 + 240 x0*x7 + 960 x0*x8 + 480 x0*x9
      + 228 x1^2 + 304 x1*x2 + 152 x1*x3 - 296 x1*x4 - 148 x1*x5 + 240 x1*x6
      + 120 x1*x7 + 480 x1*x8 + 240 x1*x9 + 448 x2^2 + 448 x2*x3 + 1072 x2*x4
      + 536 x2*x5 + 800 x2*x6 + 400 x2*x7 + 576 x2*x8 + 288 x2*x9 + 112 x3^2
      + 536 x3*x4 + 268 x3*x5 + 400 x3*x6 + 200 x3*x7 + 288 x3*x8 + 144 x3*x9
      + 4880 x4^2 + 4880 x4*x5 + 1680 x4*x6 + 840 x4*x7 + 256 x4*x8 + 128 x4*x9
      + 1220 x5^2 + 840 x5*x6 + 420 x5*x7 + 128 x5*x8 + 64 x5*x9 + 2656 x6^2
      + 2656 x6*x7 + 1120 x6*x8 + 560 x6*x9 + 664 x7^2 + 560 x7*x8 + 280 x7*x9
      + 680 x8^2 + 680 x8*x9 + 170 x9^2 ]/2
Subject To
 c0: 2 x0 + x1 + 2 x2 + x3 + 2 x4 + x5 + 2 x6 + x7 + 2 x8 + x9 = 3

Bounds
 0 <= x0 <= 1
 0 <= x1 <= 1
 0 <= x2 <= 1
 0 <= x3 <= 1
 0 <= x4 <= 1
 0 <= x5 <= 1
 0 <= x6 <= 1
 0 <= x7 <= 1
 0 <= x8 <= 1
 0 <= x9 <= 1

Binaries
 x0 x1 x2 x3 x4 x5 x6 x7 x8 x9
End
