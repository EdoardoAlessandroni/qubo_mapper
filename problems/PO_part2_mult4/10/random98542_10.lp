\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: - 200 x0 - 100 x1 + 226 x2 + 113 x3 - 402 x4 - 201 x5 + 584 x6 + 292 x7
      + 70 x8 + 35 x9 + [ 248 x0^2 + 248 x0*x1 - 224 x0*x2 - 112 x0*x3
      + 160 x0*x4 + 80 x0*x5 - 224 x0*x6 - 112 x0*x7 + 144 x0*x8 + 72 x0*x9
      + 62 x1^2 - 112 x1*x2 - 56 x1*x3 + 80 x1*x4 + 40 x1*x5 - 112 x1*x6
      - 56 x1*x7 + 72 x1*x8 + 36 x1*x9 + 832 x2^2 + 832 x2*x3 + 304 x2*x4
      + 152 x2*x5 + 928 x2*x6 + 464 x2*x7 + 160 x2*x8 + 80 x2*x9 + 208 x3^2
      + 152 x3*x4 + 76 x3*x5 + 464 x3*x6 + 232 x3*x7 + 80 x3*x8 + 40 x3*x9
      + 424 x4^2 + 424 x4*x5 - 64 x4*x6 - 32 x4*x7 + 400 x4*x8 + 200 x4*x9
      + 106 x5^2 - 32 x5*x6 - 16 x5*x7 + 200 x5*x8 + 100 x5*x9 + 1304 x6^2
      + 1304 x6*x7 + 448 x6*x8 + 224 x6*x9 + 326 x7^2 + 224 x7*x8 + 112 x7*x9
      + 440 x8^2 + 440 x8*x9 + 110 x9^2 ]/2
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
