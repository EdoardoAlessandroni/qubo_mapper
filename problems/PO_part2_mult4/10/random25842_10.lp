\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: 200 x0 + 100 x1 - 430 x2 - 215 x3 - 242 x4 - 121 x5 - 420 x6 - 210 x7
      - 78 x8 - 39 x9 + [ 808 x0^2 + 808 x0*x1 - 320 x0*x2 - 160 x0*x3
      + 912 x0*x4 + 456 x0*x5 + 496 x0*x6 + 248 x0*x7 + 576 x0*x8 + 288 x0*x9
      + 202 x1^2 - 160 x1*x2 - 80 x1*x3 + 456 x1*x4 + 228 x1*x5 + 248 x1*x6
      + 124 x1*x7 + 288 x1*x8 + 144 x1*x9 + 512 x2^2 + 512 x2*x3 + 32 x2*x4
      + 16 x2*x5 + 176 x2*x6 + 88 x2*x7 - 64 x2*x8 - 32 x2*x9 + 128 x3^2
      + 16 x3*x4 + 8 x3*x5 + 88 x3*x6 + 44 x3*x7 - 32 x3*x8 - 16 x3*x9
      + 672 x4^2 + 672 x4*x5 + 400 x4*x6 + 200 x4*x7 + 912 x4*x8 + 456 x4*x9
      + 168 x5^2 + 200 x5*x6 + 100 x5*x7 + 456 x5*x8 + 228 x5*x9 + 376 x6^2
      + 376 x6*x7 + 112 x6*x8 + 56 x6*x9 + 94 x7^2 + 56 x7*x8 + 28 x7*x9
      + 976 x8^2 + 976 x8*x9 + 244 x9^2 ]/2
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
