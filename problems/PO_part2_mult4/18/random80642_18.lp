\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: - 296 x0 - 148 x1 + 10 x2 + 5 x3 + 70 x4 + 35 x5 + 562 x6 + 281 x7
      + 208 x8 + 104 x9 + 200 x10 + 100 x11 + 226 x12 + 113 x13 - 288 x14
      - 144 x15 - 300 x16 - 150 x17 + [ 536 x0^2 + 536 x0*x1 + 112 x0*x2
      + 56 x0*x3 + 432 x0*x4 + 216 x0*x5 + 624 x0*x6 + 312 x0*x7 + 784 x0*x8
      + 392 x0*x9 + 352 x0*x10 + 176 x0*x11 + 704 x0*x12 + 352 x0*x13
      - 48 x0*x14 - 24 x0*x15 + 832 x0*x16 + 416 x0*x17 + 134 x1^2 + 56 x1*x2
      + 28 x1*x3 + 216 x1*x4 + 108 x1*x5 + 312 x1*x6 + 156 x1*x7 + 392 x1*x8
      + 196 x1*x9 + 176 x1*x10 + 88 x1*x11 + 352 x1*x12 + 176 x1*x13 - 24 x1*x14
      - 12 x1*x15 + 416 x1*x16 + 208 x1*x17 + 1000 x2^2 + 1000 x2*x3 + 64 x2*x4
      + 32 x2*x5 + 16 x2*x6 + 8 x2*x7 + 128 x2*x8 + 64 x2*x9 + 1072 x2*x10
      + 536 x2*x11 - 64 x2*x12 - 32 x2*x13 - 48 x2*x14 - 24 x2*x15 + 128 x2*x16
      + 64 x2*x17 + 250 x3^2 + 32 x3*x4 + 16 x3*x5 + 8 x3*x6 + 4 x3*x7
      + 64 x3*x8 + 32 x3*x9 + 536 x3*x10 + 268 x3*x11 - 32 x3*x12 - 16 x3*x13
      - 24 x3*x14 - 12 x3*x15 + 64 x3*x16 + 32 x3*x17 + 440 x4^2 + 440 x4*x5
      + 320 x4*x6 + 160 x4*x7 + 240 x4*x8 + 120 x4*x9 + 160 x4*x10 + 80 x4*x11
      + 160 x4*x12 + 80 x4*x13 + 272 x4*x14 + 136 x4*x15 + 624 x4*x16
      + 312 x4*x17 + 110 x5^2 + 160 x5*x6 + 80 x5*x7 + 120 x5*x8 + 60 x5*x9
      + 80 x5*x10 + 40 x5*x11 + 80 x5*x12 + 40 x5*x13 + 136 x5*x14 + 68 x5*x15
      + 312 x5*x16 + 156 x5*x17 + 608 x6^2 + 608 x6*x7 + 592 x6*x8 + 296 x6*x9
      + 48 x6*x10 + 24 x6*x11 + 320 x6*x12 + 160 x6*x13 + 64 x6*x14 + 32 x6*x15
      + 32 x6*x16 + 16 x6*x17 + 152 x7^2 + 296 x7*x8 + 148 x7*x9 + 24 x7*x10
      + 12 x7*x11 + 160 x7*x12 + 80 x7*x13 + 32 x7*x14 + 16 x7*x15 + 16 x7*x16
      + 8 x7*x17 + 680 x8^2 + 680 x8*x9 + 352 x8*x10 + 176 x8*x11 + 896 x8*x12
      + 448 x8*x13 - 64 x8*x14 - 32 x8*x15 + 800 x8*x16 + 400 x8*x17 + 170 x9^2
      + 176 x9*x10 + 88 x9*x11 + 448 x9*x12 + 224 x9*x13 - 32 x9*x14 - 16 x9*x15
      + 400 x9*x16 + 200 x9*x17 + 808 x10^2 + 808 x10*x11 + 688 x10*x12
      + 344 x10*x13 - 80 x10*x14 - 40 x10*x15 + 464 x10*x16 + 232 x10*x17
      + 202 x11^2 + 344 x11*x12 + 172 x11*x13 - 40 x11*x14 - 20 x11*x15
      + 232 x11*x16 + 116 x11*x17 + 832 x12^2 + 832 x12*x13 - 384 x12*x14
      - 192 x12*x15 + 1440 x12*x16 + 720 x12*x17 + 208 x13^2 - 192 x13*x14
      - 96 x13*x15 + 720 x13*x16 + 360 x13*x17 + 320 x14^2 + 320 x14*x15
      + 112 x14*x16 + 56 x14*x17 + 80 x15^2 + 56 x15*x16 + 28 x15*x17
      + 1600 x16^2 + 1600 x16*x17 + 400 x17^2 ]/2
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
