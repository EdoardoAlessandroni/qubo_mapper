\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: 200 x0 + 100 x1 + 56 x2 + 28 x3 + 90 x4 + 45 x5 - 128 x6 - 64 x7 - 88 x8
      - 44 x9 + [ 808 x0^2 + 808 x0*x1 + 160 x0*x2 + 80 x0*x3 + 464 x0*x4
      + 232 x0*x5 + 80 x0*x6 + 40 x0*x7 + 448 x0*x8 + 224 x0*x9 + 202 x1^2
      + 80 x1*x2 + 40 x1*x3 + 232 x1*x4 + 116 x1*x5 + 40 x1*x6 + 20 x1*x7
      + 224 x1*x8 + 112 x1*x9 + 536 x2^2 + 536 x2*x3 + 288 x2*x4 + 144 x2*x5
      + 320 x2*x6 + 160 x2*x7 + 560 x2*x8 + 280 x2*x9 + 134 x3^2 + 144 x3*x4
      + 72 x3*x5 + 160 x3*x6 + 80 x3*x7 + 280 x3*x8 + 140 x3*x9 + 664 x4^2
      + 664 x4*x5 + 272 x4*x6 + 136 x4*x7 + 544 x4*x8 + 272 x4*x9 + 166 x5^2
      + 136 x5*x6 + 68 x5*x7 + 272 x5*x8 + 136 x5*x9 + 368 x6^2 + 368 x6*x7
      + 384 x6*x8 + 192 x6*x9 + 92 x7^2 + 192 x7*x8 + 96 x7*x9 + 472 x8^2
      + 472 x8*x9 + 118 x9^2 ]/2
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
