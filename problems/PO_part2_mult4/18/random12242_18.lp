\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: 308 x0 + 154 x1 - 224 x2 - 112 x3 + 206 x4 + 103 x5 - 200 x6 - 100 x7
      - 430 x8 - 215 x9 - 52 x10 - 26 x11 - 152 x12 - 76 x13 - 288 x14 - 144 x15
      - 672 x16 - 336 x17 + [ 1128 x0^2 + 1128 x0*x1 + 608 x0*x2 + 304 x0*x3
      + 1248 x0*x4 + 624 x0*x5 + 192 x0*x6 + 96 x0*x7 + 544 x0*x8 + 272 x0*x9
      + 512 x0*x10 + 256 x0*x11 + 720 x0*x12 + 360 x0*x13 + 64 x0*x14
      + 32 x0*x15 + 1184 x0*x16 + 592 x0*x17 + 282 x1^2 + 304 x1*x2 + 152 x1*x3
      + 624 x1*x4 + 312 x1*x5 + 96 x1*x6 + 48 x1*x7 + 272 x1*x8 + 136 x1*x9
      + 256 x1*x10 + 128 x1*x11 + 360 x1*x12 + 180 x1*x13 + 32 x1*x14
      + 16 x1*x15 + 592 x1*x16 + 296 x1*x17 + 1152 x2^2 + 1152 x2*x3 + 976 x2*x4
      + 488 x2*x5 + 320 x2*x6 + 160 x2*x7 + 448 x2*x8 + 224 x2*x9 + 544 x2*x10
      + 272 x2*x11 + 400 x2*x12 + 200 x2*x13 - 48 x2*x14 - 24 x2*x15
      + 464 x2*x16 + 232 x2*x17 + 288 x3^2 + 488 x3*x4 + 244 x3*x5 + 160 x3*x6
      + 80 x3*x7 + 224 x3*x8 + 112 x3*x9 + 272 x3*x10 + 136 x3*x11 + 200 x3*x12
      + 100 x3*x13 - 24 x3*x14 - 12 x3*x15 + 232 x3*x16 + 116 x3*x17 + 984 x4^2
      + 984 x4*x5 + 368 x4*x6 + 184 x4*x7 + 272 x4*x8 + 136 x4*x9 + 384 x4*x10
      + 192 x4*x11 + 768 x4*x12 + 384 x4*x13 + 128 x4*x14 + 64 x4*x15
      + 896 x4*x16 + 448 x4*x17 + 246 x5^2 + 184 x5*x6 + 92 x5*x7 + 136 x5*x8
      + 68 x5*x9 + 192 x5*x10 + 96 x5*x11 + 384 x5*x12 + 192 x5*x13 + 64 x5*x14
      + 32 x5*x15 + 448 x5*x16 + 224 x5*x17 + 248 x6^2 + 248 x6*x7 + 80 x6*x8
      + 40 x6*x9 - 96 x6*x10 - 48 x6*x11 + 80 x6*x12 + 40 x6*x13 + 304 x6*x14
      + 152 x6*x15 + 32 x6*x16 + 16 x6*x17 + 62 x7^2 + 40 x7*x8 + 20 x7*x9
      - 48 x7*x10 - 24 x7*x11 + 40 x7*x12 + 20 x7*x13 + 152 x7*x14 + 76 x7*x15
      + 16 x7*x16 + 8 x7*x17 + 512 x8^2 + 512 x8*x9 + 416 x8*x10 + 208 x8*x11
      + 144 x8*x12 + 72 x8*x13 + 144 x8*x14 + 72 x8*x15 + 384 x8*x16
      + 192 x8*x17 + 128 x9^2 + 208 x9*x10 + 104 x9*x11 + 72 x9*x12 + 36 x9*x13
      + 72 x9*x14 + 36 x9*x15 + 192 x9*x16 + 96 x9*x17 + 376 x10^2 + 376 x10*x11
      + 256 x10*x12 + 128 x10*x13 + 64 x10*x14 + 32 x10*x15 + 592 x10*x16
      + 296 x10*x17 + 94 x11^2 + 128 x11*x12 + 64 x11*x13 + 32 x11*x14
      + 16 x11*x15 + 296 x11*x16 + 148 x11*x17 + 480 x12^2 + 480 x12*x13
      - 112 x12*x14 - 56 x12*x15 + 672 x12*x16 + 336 x12*x17 + 120 x13^2
      - 56 x13*x14 - 28 x13*x15 + 336 x13*x16 + 168 x13*x17 + 320 x14^2
      + 320 x14*x15 + 208 x14*x16 + 104 x14*x17 + 80 x15^2 + 104 x15*x16
      + 52 x15*x17 + 1048 x16^2 + 1048 x16*x17 + 262 x17^2 ]/2
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
