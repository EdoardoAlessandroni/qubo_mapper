\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: 364 x0 + 182 x1 - 590 x2 - 295 x3 + 230 x4 + 115 x5 - 420 x6 - 210 x7
      + 584 x8 + 292 x9 + [ 608 x0^2 + 608 x0*x1 + 144 x0*x2 + 72 x0*x3
      - 96 x0*x4 - 48 x0*x5 + 256 x0*x6 + 128 x0*x7 + 704 x0*x8 + 352 x0*x9
      + 152 x1^2 + 72 x1*x2 + 36 x1*x3 - 48 x1*x4 - 24 x1*x5 + 128 x1*x6
      + 64 x1*x7 + 352 x1*x8 + 176 x1*x9 + 392 x2^2 + 392 x2*x3 + 608 x2*x4
      + 304 x2*x5 + 160 x2*x6 + 80 x2*x7 + 96 x2*x8 + 48 x2*x9 + 98 x3^2
      + 304 x3*x4 + 152 x3*x5 + 80 x3*x6 + 40 x3*x7 + 48 x3*x8 + 24 x3*x9
      + 1000 x4^2 + 1000 x4*x5 + 384 x4*x6 + 192 x4*x7 + 720 x4*x8 + 360 x4*x9
      + 250 x5^2 + 192 x5*x6 + 96 x5*x7 + 360 x5*x8 + 180 x5*x9 + 376 x6^2
      + 376 x6*x7 + 400 x6*x8 + 200 x6*x9 + 94 x7^2 + 200 x7*x8 + 100 x7*x9
      + 1304 x8^2 + 1304 x8*x9 + 326 x9^2 ]/2
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
