\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: 374 x0 + 187 x1 + 586 x2 + 293 x3 - 590 x4 - 295 x5 - 106 x6 - 53 x7
      - 108 x8 - 54 x9 + 584 x10 + 292 x11 + 686 x12 + 343 x13 + 364 x14
      + 182 x15 - 154 x16 - 77 x17 + [ 600 x0^2 + 600 x0*x1 + 96 x0*x2
      + 48 x0*x3 + 496 x0*x4 + 248 x0*x5 + 496 x0*x6 + 248 x0*x7 + 656 x0*x8
      + 328 x0*x9 + 512 x0*x10 + 256 x0*x11 + 848 x0*x12 + 424 x0*x13
      + 448 x0*x14 + 224 x0*x15 + 704 x0*x16 + 352 x0*x17 + 150 x1^2 + 48 x1*x2
      + 24 x1*x3 + 248 x1*x4 + 124 x1*x5 + 248 x1*x6 + 124 x1*x7 + 328 x1*x8
      + 164 x1*x9 + 256 x1*x10 + 128 x1*x11 + 424 x1*x12 + 212 x1*x13
      + 224 x1*x14 + 112 x1*x15 + 352 x1*x16 + 176 x1*x17 + 1240 x2^2
      + 1240 x2*x3 + 384 x2*x4 + 192 x2*x5 + 96 x2*x6 + 48 x2*x7 + 256 x2*x8
      + 128 x2*x9 - 80 x2*x10 - 40 x2*x11 + 448 x2*x12 + 224 x2*x13 + 128 x2*x14
      + 64 x2*x15 - 512 x2*x16 - 256 x2*x17 + 310 x3^2 + 192 x3*x4 + 96 x3*x5
      + 48 x3*x6 + 24 x3*x7 + 128 x3*x8 + 64 x3*x9 - 40 x3*x10 - 20 x3*x11
      + 224 x3*x12 + 112 x3*x13 + 64 x3*x14 + 32 x3*x15 - 256 x3*x16
      - 128 x3*x17 + 392 x4^2 + 392 x4*x5 + 336 x4*x6 + 168 x4*x7 + 384 x4*x8
      + 192 x4*x9 + 96 x4*x10 + 48 x4*x11 + 576 x4*x12 + 288 x4*x13 + 144 x4*x14
      + 72 x4*x15 + 496 x4*x16 + 248 x4*x17 + 98 x5^2 + 168 x5*x6 + 84 x5*x7
      + 192 x5*x8 + 96 x5*x9 + 48 x5*x10 + 24 x5*x11 + 288 x5*x12 + 144 x5*x13
      + 72 x5*x14 + 36 x5*x15 + 248 x5*x16 + 124 x5*x17 + 344 x6^2 + 344 x6*x7
      + 432 x6*x8 + 216 x6*x9 + 160 x6*x10 + 80 x6*x11 + 176 x6*x12 + 88 x6*x13
      + 272 x6*x14 + 136 x6*x15 + 272 x6*x16 + 136 x6*x17 + 86 x7^2 + 216 x7*x8
      + 108 x7*x9 + 80 x7*x10 + 40 x7*x11 + 88 x7*x12 + 44 x7*x13 + 136 x7*x14
      + 68 x7*x15 + 136 x7*x16 + 68 x7*x17 + 360 x8^2 + 360 x8*x9 + 432 x8*x10
      + 216 x8*x11 + 528 x8*x12 + 264 x8*x13 + 352 x8*x14 + 176 x8*x15
      + 240 x8*x16 + 120 x8*x17 + 90 x9^2 + 216 x9*x10 + 108 x9*x11 + 264 x9*x12
      + 132 x9*x13 + 176 x9*x14 + 88 x9*x15 + 120 x9*x16 + 60 x9*x17
      + 1304 x10^2 + 1304 x10*x11 + 1456 x10*x12 + 728 x10*x13 + 704 x10*x14
      + 352 x10*x15 - 384 x10*x16 - 192 x10*x17 + 326 x11^2 + 728 x11*x12
      + 364 x11*x13 + 352 x11*x14 + 176 x11*x15 - 192 x11*x16 - 96 x11*x17
      + 1352 x12^2 + 1352 x12*x13 + 560 x12*x14 + 280 x12*x15 + 448 x12*x16
      + 224 x12*x17 + 338 x13^2 + 280 x13*x14 + 140 x13*x15 + 224 x13*x16
      + 112 x13*x17 + 608 x14^2 + 608 x14*x15 - 32 x14*x16 - 16 x14*x17
      + 152 x15^2 - 16 x15*x16 - 8 x15*x17 + 1264 x16^2 + 1264 x16*x17
      + 316 x17^2 ]/2
Subject To
 c0: 2 x0 + x1 + 2 x2 + x3 + 2 x4 + x5 + 2 x6 + x7 + 2 x8 + x9 + 2 x10 + x11
     + 2 x12 + x13 + 2 x14 + x15 + 2 x16 + x17 = 3

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
 0 <= x10 <= 1
 0 <= x11 <= 1
 0 <= x12 <= 1
 0 <= x13 <= 1
 0 <= x14 <= 1
 0 <= x15 <= 1
 0 <= x16 <= 1
 0 <= x17 <= 1

Binaries
 x0 x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17
End
