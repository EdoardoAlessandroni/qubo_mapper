\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: - 672 x0 - 336 x1 + 374 x2 + 187 x3 + 370 x4 + 185 x5 + 56 x6 + 28 x7
      - 422 x8 - 211 x9 - 86 x10 - 43 x11 + 110 x12 + 55 x13 + 252 x14 + 126 x15
      - 156 x16 - 78 x17 + [ 1048 x0^2 + 1048 x0*x1 + 720 x0*x2 + 360 x0*x3
      + 272 x0*x4 + 136 x0*x5 + 688 x0*x6 + 344 x0*x7 + 448 x0*x8 + 224 x0*x9
      + 448 x0*x10 + 224 x0*x11 + 656 x0*x12 + 328 x0*x13 - 64 x0*x14
      - 32 x0*x15 + 672 x0*x16 + 336 x0*x17 + 262 x1^2 + 360 x1*x2 + 180 x1*x3
      + 136 x1*x4 + 68 x1*x5 + 344 x1*x6 + 172 x1*x7 + 224 x1*x8 + 112 x1*x9
      + 224 x1*x10 + 112 x1*x11 + 328 x1*x12 + 164 x1*x13 - 32 x1*x14
      - 16 x1*x15 + 336 x1*x16 + 168 x1*x17 + 600 x2^2 + 600 x2*x3 + 768 x2*x4
      + 384 x2*x5 + 400 x2*x6 + 200 x2*x7 + 880 x2*x8 + 440 x2*x9 + 832 x2*x10
      + 416 x2*x11 + 624 x2*x12 + 312 x2*x13 + 480 x2*x14 + 240 x2*x15
      + 640 x2*x16 + 320 x2*x17 + 150 x3^2 + 384 x3*x4 + 192 x3*x5 + 200 x3*x6
      + 100 x3*x7 + 440 x3*x8 + 220 x3*x9 + 416 x3*x10 + 208 x3*x11 + 312 x3*x12
      + 156 x3*x13 + 240 x3*x14 + 120 x3*x15 + 320 x3*x16 + 160 x3*x17
      + 960 x4^2 + 960 x4*x5 + 224 x4*x6 + 112 x4*x7 + 976 x4*x8 + 488 x4*x9
      + 320 x4*x10 + 160 x4*x11 - 48 x4*x12 - 24 x4*x13 + 560 x4*x14
      + 280 x4*x15 + 720 x4*x16 + 360 x4*x17 + 240 x5^2 + 112 x5*x6 + 56 x5*x7
      + 488 x5*x8 + 244 x5*x9 + 160 x5*x10 + 80 x5*x11 - 24 x5*x12 - 12 x5*x13
      + 280 x5*x14 + 140 x5*x15 + 360 x5*x16 + 180 x5*x17 + 536 x6^2 + 536 x6*x7
      + 384 x6*x8 + 192 x6*x9 + 480 x6*x10 + 240 x6*x11 + 144 x6*x12 + 72 x6*x13
      + 256 x6*x14 + 128 x6*x15 + 240 x6*x16 + 120 x6*x17 + 134 x7^2 + 192 x7*x8
      + 96 x7*x9 + 240 x7*x10 + 120 x7*x11 + 72 x7*x12 + 36 x7*x13 + 128 x7*x14
      + 64 x7*x15 + 120 x7*x16 + 60 x7*x17 + 576 x8^2 + 576 x8*x9 + 704 x8*x10
      + 352 x8*x11 + 544 x8*x12 + 272 x8*x13 + 592 x8*x14 + 296 x8*x15
      + 720 x8*x16 + 360 x8*x17 + 144 x9^2 + 352 x9*x10 + 176 x9*x11
      + 272 x9*x12 + 136 x9*x13 + 296 x9*x14 + 148 x9*x15 + 360 x9*x16
      + 180 x9*x17 + 496 x10^2 + 496 x10*x11 + 320 x10*x12 + 160 x10*x13
      + 400 x10*x14 + 200 x10*x15 + 272 x10*x16 + 136 x10*x17 + 124 x11^2
      + 160 x11*x12 + 80 x11*x13 + 200 x11*x14 + 100 x11*x15 + 136 x11*x16
      + 68 x11*x17 + 1192 x12^2 + 1192 x12*x13 + 912 x12*x14 + 456 x12*x15
      + 496 x12*x16 + 248 x12*x17 + 298 x13^2 + 456 x13*x14 + 228 x13*x15
      + 248 x13*x16 + 124 x13*x17 + 536 x14^2 + 536 x14*x15 + 400 x14*x16
      + 200 x14*x17 + 134 x15^2 + 200 x15*x16 + 100 x15*x17 + 504 x16^2
      + 504 x16*x17 + 126 x17^2 ]/2
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
