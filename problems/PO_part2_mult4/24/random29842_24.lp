\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: - 196 x0 - 98 x1 + 90 x2 + 45 x3 + 498 x4 + 249 x5 + 206 x6 + 103 x7
      - 110 x8 - 55 x9 + 54 x10 + 27 x11 - 10 x12 - 5 x13 - 228 x14 - 114 x15
      + 400 x16 + 200 x17 + 376 x18 + 188 x19 + 200 x20 + 100 x21 + 72 x22
      + 36 x23 + [ 344 x0^2 + 344 x0*x1 + 192 x0*x2 + 96 x0*x3 + 752 x0*x4
      + 376 x0*x5 + 848 x0*x6 + 424 x0*x7 + 944 x0*x8 + 472 x0*x9 + 160 x0*x10
      + 80 x0*x11 + 720 x0*x12 + 360 x0*x13 + 608 x0*x14 + 304 x0*x15
      + 1264 x0*x16 + 632 x0*x17 + 384 x0*x18 + 192 x0*x19 + 496 x0*x20
      + 248 x0*x21 + 112 x0*x22 + 56 x0*x23 + 86 x1^2 + 96 x1*x2 + 48 x1*x3
      + 376 x1*x4 + 188 x1*x5 + 424 x1*x6 + 212 x1*x7 + 472 x1*x8 + 236 x1*x9
      + 80 x1*x10 + 40 x1*x11 + 360 x1*x12 + 180 x1*x13 + 304 x1*x14
      + 152 x1*x15 + 632 x1*x16 + 316 x1*x17 + 192 x1*x18 + 96 x1*x19
      + 248 x1*x20 + 124 x1*x21 + 56 x1*x22 + 28 x1*x23 + 2528 x2^2 + 2528 x2*x3
      + 784 x2*x4 + 392 x2*x5 + 976 x2*x6 + 488 x2*x7 + 1344 x2*x8 + 672 x2*x9
      - 208 x2*x10 - 104 x2*x11 + 192 x2*x12 + 96 x2*x13 + 800 x2*x14
      + 400 x2*x15 + 624 x2*x16 + 312 x2*x17 + 1808 x2*x18 + 904 x2*x19
      + 304 x2*x20 + 152 x2*x21 + 1056 x2*x22 + 528 x2*x23 + 632 x3^2
      + 392 x3*x4 + 196 x3*x5 + 488 x3*x6 + 244 x3*x7 + 672 x3*x8 + 336 x3*x9
      - 104 x3*x10 - 52 x3*x11 + 96 x3*x12 + 48 x3*x13 + 400 x3*x14 + 200 x3*x15
      + 312 x3*x16 + 156 x3*x17 + 904 x3*x18 + 452 x3*x19 + 152 x3*x20
      + 76 x3*x21 + 528 x3*x22 + 264 x3*x23 + 680 x4^2 + 680 x4*x5 + 1120 x4*x6
      + 560 x4*x7 + 1200 x4*x8 + 600 x4*x9 + 960 x4*x12 + 480 x4*x13
      + 704 x4*x14 + 352 x4*x15 + 1888 x4*x16 + 944 x4*x17 + 704 x4*x18
      + 352 x4*x19 + 496 x4*x20 + 248 x4*x21 + 320 x4*x22 + 160 x4*x23
      + 170 x5^2 + 560 x5*x6 + 280 x5*x7 + 600 x5*x8 + 300 x5*x9 + 480 x5*x12
      + 240 x5*x13 + 352 x5*x14 + 176 x5*x15 + 944 x5*x16 + 472 x5*x17
      + 352 x5*x18 + 176 x5*x19 + 248 x5*x20 + 124 x5*x21 + 160 x5*x22
      + 80 x5*x23 + 984 x6^2 + 984 x6*x7 + 1664 x6*x8 + 832 x6*x9 + 352 x6*x10
      + 176 x6*x11 + 576 x6*x12 + 288 x6*x13 + 1008 x6*x14 + 504 x6*x15
      + 1120 x6*x16 + 560 x6*x17 + 848 x6*x18 + 424 x6*x19 + 704 x6*x20
      + 352 x6*x21 + 576 x6*x22 + 288 x6*x23 + 246 x7^2 + 832 x7*x8 + 416 x7*x9
      + 176 x7*x10 + 88 x7*x11 + 288 x7*x12 + 144 x7*x13 + 504 x7*x14
      + 252 x7*x15 + 560 x7*x16 + 280 x7*x17 + 424 x7*x18 + 212 x7*x19
      + 352 x7*x20 + 176 x7*x21 + 288 x7*x22 + 144 x7*x23 + 1224 x8^2
      + 1224 x8*x9 + 608 x8*x10 + 304 x8*x11 + 672 x8*x12 + 336 x8*x13
      + 1408 x8*x14 + 704 x8*x15 + 1664 x8*x16 + 832 x8*x17 + 1008 x8*x18
      + 504 x8*x19 + 592 x8*x20 + 296 x8*x21 + 704 x8*x22 + 352 x8*x23
      + 306 x9^2 + 304 x9*x10 + 152 x9*x11 + 336 x9*x12 + 168 x9*x13
      + 704 x9*x14 + 352 x9*x15 + 832 x9*x16 + 416 x9*x17 + 504 x9*x18
      + 252 x9*x19 + 296 x9*x20 + 148 x9*x21 + 352 x9*x22 + 176 x9*x23
      + 392 x10^2 + 392 x10*x11 + 416 x10*x14 + 208 x10*x15 + 80 x10*x16
      + 40 x10*x17 + 272 x10*x18 + 136 x10*x19 + 32 x10*x20 + 16 x10*x21
      + 256 x10*x22 + 128 x10*x23 + 98 x11^2 + 208 x11*x14 + 104 x11*x15
      + 40 x11*x16 + 20 x11*x17 + 136 x11*x18 + 68 x11*x19 + 16 x11*x20
      + 8 x11*x21 + 128 x11*x22 + 64 x11*x23 + 912 x12^2 + 912 x12*x13
      + 496 x12*x14 + 248 x12*x15 + 2736 x12*x16 + 1368 x12*x17 + 224 x12*x18
      + 112 x12*x19 + 928 x12*x20 + 464 x12*x21 + 112 x12*x22 + 56 x12*x23
      + 228 x13^2 + 248 x13*x14 + 124 x13*x15 + 1368 x13*x16 + 684 x13*x17
      + 112 x13*x18 + 56 x13*x19 + 464 x13*x20 + 232 x13*x21 + 56 x13*x22
      + 28 x13*x23 + 632 x14^2 + 632 x14*x15 + 1088 x14*x16 + 544 x14*x17
      + 672 x14*x18 + 336 x14*x19 + 880 x14*x20 + 440 x14*x21 + 368 x14*x22
      + 184 x14*x23 + 158 x15^2 + 544 x15*x16 + 272 x15*x17 + 336 x15*x18
      + 168 x15*x19 + 440 x15*x20 + 220 x15*x21 + 184 x15*x22 + 92 x15*x23
      + 2528 x16^2 + 2528 x16*x17 + 416 x16*x18 + 208 x16*x19 + 1456 x16*x20
      + 728 x16*x21 + 48 x16*x22 + 24 x16*x23 + 632 x17^2 + 208 x17*x18
      + 104 x17*x19 + 728 x17*x20 + 364 x17*x21 + 24 x17*x22 + 12 x17*x23
      + 792 x18^2 + 792 x18*x19 + 16 x18*x20 + 8 x18*x21 + 672 x18*x22
      + 336 x18*x23 + 198 x19^2 + 8 x19*x20 + 4 x19*x21 + 336 x19*x22
      + 168 x19*x23 + 808 x20^2 + 808 x20*x21 + 208 x20*x22 + 104 x20*x23
      + 202 x21^2 + 104 x21*x22 + 52 x21*x23 + 552 x22^2 + 552 x22*x23
      + 138 x23^2 ]/2
Subject To
 c0: 2 x0 + x1 + 2 x2 + x3 + 2 x4 + x5 + 2 x6 + x7 + 2 x8 + x9 + 2 x10 + x11
     + 2 x12 + x13 + 2 x14 + x15 + 2 x16 + x17 + 2 x18 + x19 + 2 x20 + x21
     + 2 x22 + x23 = 3

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
 0 <= x18 <= 1
 0 <= x19 <= 1
 0 <= x20 <= 1
 0 <= x21 <= 1
 0 <= x22 <= 1
 0 <= x23 <= 1

Binaries
 x0 x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18 x19 x20 x21
 x22 x23
End
