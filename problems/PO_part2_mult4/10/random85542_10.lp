\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: 230 x0 + 115 x1 - 52 x2 - 26 x3 + 70 x4 + 35 x5 + 42 x6 + 21 x7 - 340 x8
      - 170 x9 + [ 1000 x0^2 + 1000 x0*x1 + 464 x0*x2 + 232 x0*x3 + 640 x0*x4
      + 320 x0*x5 + 896 x0*x6 + 448 x0*x7 + 352 x0*x8 + 176 x0*x9 + 250 x1^2
      + 232 x1*x2 + 116 x1*x3 + 320 x1*x4 + 160 x1*x5 + 448 x1*x6 + 224 x1*x7
      + 176 x1*x8 + 88 x1*x9 + 376 x2^2 + 376 x2*x3 + 400 x2*x4 + 200 x2*x5
      + 224 x2*x6 + 112 x2*x7 + 96 x2*x8 + 48 x2*x9 + 94 x3^2 + 200 x3*x4
      + 100 x3*x5 + 112 x3*x6 + 56 x3*x7 + 48 x3*x8 + 24 x3*x9 + 440 x4^2
      + 440 x4*x5 + 464 x4*x6 + 232 x4*x7 + 432 x4*x8 + 216 x4*x9 + 110 x5^2
      + 232 x5*x6 + 116 x5*x7 + 216 x5*x8 + 108 x5*x9 + 1144 x6^2 + 1144 x6*x7
      + 720 x6*x8 + 360 x6*x9 + 286 x7^2 + 360 x7*x8 + 180 x7*x9 + 528 x8^2
      + 528 x8*x9 + 132 x9^2 ]/2
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
