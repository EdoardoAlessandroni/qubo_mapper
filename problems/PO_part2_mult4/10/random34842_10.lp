\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: 40 x0 + 20 x1 + 200 x2 + 100 x3 - 68 x4 - 34 x5 + 110 x6 + 55 x7 - 672 x8
      - 336 x9 + [ 984 x0^2 + 984 x0*x1 + 592 x0*x2 + 296 x0*x3 + 1584 x0*x4
      + 792 x0*x5 + 912 x0*x6 + 456 x0*x7 + 1216 x0*x8 + 608 x0*x9 + 246 x1^2
      + 296 x1*x2 + 148 x1*x3 + 792 x1*x4 + 396 x1*x5 + 456 x1*x6 + 228 x1*x7
      + 608 x1*x8 + 304 x1*x9 + 808 x2^2 + 808 x2*x3 + 528 x2*x4 + 264 x2*x5
      + 784 x2*x6 + 392 x2*x7 + 144 x2*x8 + 72 x2*x9 + 202 x3^2 + 264 x3*x4
      + 132 x3*x5 + 392 x3*x6 + 196 x3*x7 + 72 x3*x8 + 36 x3*x9 + 2160 x4^2
      + 2160 x4*x5 + 432 x4*x6 + 216 x4*x7 + 1728 x4*x8 + 864 x4*x9 + 540 x5^2
      + 216 x5*x6 + 108 x5*x7 + 864 x5*x8 + 432 x5*x9 + 1192 x6^2 + 1192 x6*x7
      + 656 x6*x8 + 328 x6*x9 + 298 x7^2 + 328 x7*x8 + 164 x7*x9 + 1048 x8^2
      + 1048 x8*x9 + 262 x9^2 ]/2
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
