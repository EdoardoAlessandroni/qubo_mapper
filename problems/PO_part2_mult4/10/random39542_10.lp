\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: 258 x0 + 129 x1 - 1290 x2 - 645 x3 + 52 x4 + 26 x5 + 90 x6 + 45 x7
      + 216 x8 + 108 x9 + [ 1944 x0^2 + 1944 x0*x1 + 256 x0*x2 + 128 x0*x3
      + 2176 x0*x4 + 1088 x0*x5 + 768 x0*x6 + 384 x0*x7 + 1408 x0*x8 + 704 x0*x9
      + 486 x1^2 + 128 x1*x2 + 64 x1*x3 + 1088 x1*x4 + 544 x1*x5 + 384 x1*x6
      + 192 x1*x7 + 704 x1*x8 + 352 x1*x9 + 2080 x2^2 + 2080 x2*x3 + 992 x2*x4
      + 496 x2*x5 + 256 x2*x6 + 128 x2*x7 + 624 x2*x8 + 312 x2*x9 + 520 x3^2
      + 496 x3*x4 + 248 x3*x5 + 128 x3*x6 + 64 x3*x7 + 312 x3*x8 + 156 x3*x9
      + 2448 x4^2 + 2448 x4*x5 + 1200 x4*x6 + 600 x4*x7 + 1744 x4*x8 + 872 x4*x9
      + 612 x5^2 + 600 x5*x6 + 300 x5*x7 + 872 x5*x8 + 436 x5*x9 + 448 x6^2
      + 448 x6*x7 + 640 x6*x8 + 320 x6*x9 + 112 x7^2 + 320 x7*x8 + 160 x7*x9
      + 1144 x8^2 + 1144 x8*x9 + 286 x9^2 ]/2
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
