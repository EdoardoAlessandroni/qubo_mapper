\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: 284 x0 + 142 x1 + 42 x2 + 21 x3 - 128 x4 - 64 x5 - 2 x6 - x7 - 312 x8
      - 156 x9 + 54 x10 + 27 x11 + 628 x12 + 314 x13 - 78 x14 - 39 x15 - 856 x16
      - 428 x17 - 236 x18 - 118 x19 + [ 1104 x0^2 + 1104 x0*x1 + 1280 x0*x2
      + 640 x0*x3 + 400 x0*x4 + 200 x0*x5 + 1216 x0*x6 + 608 x0*x7 + 128 x0*x8
      + 64 x0*x9 + 80 x0*x10 + 40 x0*x11 + 96 x0*x12 + 48 x0*x13 + 720 x0*x14
      + 360 x0*x15 + 1728 x0*x16 + 864 x0*x17 + 1344 x0*x18 + 672 x0*x19
      + 276 x1^2 + 640 x1*x2 + 320 x1*x3 + 200 x1*x4 + 100 x1*x5 + 608 x1*x6
      + 304 x1*x7 + 64 x1*x8 + 32 x1*x9 + 40 x1*x10 + 20 x1*x11 + 48 x1*x12
      + 24 x1*x13 + 360 x1*x14 + 180 x1*x15 + 864 x1*x16 + 432 x1*x17
      + 672 x1*x18 + 336 x1*x19 + 1016 x2^2 + 1016 x2*x3 + 256 x2*x4 + 128 x2*x5
      + 736 x2*x6 + 368 x2*x7 + 16 x2*x8 + 8 x2*x9 + 192 x2*x10 + 96 x2*x11
      + 912 x2*x12 + 456 x2*x13 + 928 x2*x14 + 464 x2*x15 + 1696 x2*x16
      + 848 x2*x17 + 1808 x2*x18 + 904 x2*x19 + 254 x3^2 + 128 x3*x4 + 64 x3*x5
      + 368 x3*x6 + 184 x3*x7 + 8 x3*x8 + 4 x3*x9 + 96 x3*x10 + 48 x3*x11
      + 456 x3*x12 + 228 x3*x13 + 464 x3*x14 + 232 x3*x15 + 848 x3*x16
      + 424 x3*x17 + 904 x3*x18 + 452 x3*x19 + 368 x4^2 + 368 x4*x5 - 368 x4*x6
      - 184 x4*x7 - 32 x4*x8 - 16 x4*x9 + 288 x4*x10 + 144 x4*x11 - 96 x4*x12
      - 48 x4*x13 + 496 x4*x14 + 248 x4*x15 + 320 x4*x16 + 160 x4*x17
      + 400 x4*x18 + 200 x4*x19 + 92 x5^2 - 184 x5*x6 - 92 x5*x7 - 16 x5*x8
      - 8 x5*x9 + 144 x5*x10 + 72 x5*x11 - 48 x5*x12 - 24 x5*x13 + 248 x5*x14
      + 124 x5*x15 + 160 x5*x16 + 80 x5*x17 + 200 x5*x18 + 100 x5*x19
      + 2128 x6^2 + 2128 x6*x7 + 96 x6*x8 + 48 x6*x9 + 336 x6*x10 + 168 x6*x11
      + 528 x6*x12 + 264 x6*x13 + 80 x6*x14 + 40 x6*x15 + 3472 x6*x16
      + 1736 x6*x17 + 1664 x6*x18 + 832 x6*x19 + 532 x7^2 + 48 x7*x8 + 24 x7*x9
      + 168 x7*x10 + 84 x7*x11 + 264 x7*x12 + 132 x7*x13 + 40 x7*x14 + 20 x7*x15
      + 1736 x7*x16 + 868 x7*x17 + 832 x7*x18 + 416 x7*x19 + 1168 x8^2
      + 1168 x8*x9 + 272 x8*x10 + 136 x8*x11 + 208 x8*x12 + 104 x8*x13
      - 256 x8*x14 - 128 x8*x15 + 512 x8*x16 + 256 x8*x17 + 292 x9^2
      + 136 x9*x10 + 68 x9*x11 + 104 x9*x12 + 52 x9*x13 - 128 x9*x14 - 64 x9*x15
      + 256 x9*x16 + 128 x9*x17 + 392 x10^2 + 392 x10*x11 + 96 x10*x12
      + 48 x10*x13 + 144 x10*x14 + 72 x10*x15 + 208 x10*x16 + 104 x10*x17
      + 112 x10*x18 + 56 x10*x19 + 98 x11^2 + 48 x11*x12 + 24 x11*x13
      + 72 x11*x14 + 36 x11*x15 + 104 x11*x16 + 52 x11*x17 + 56 x11*x18
      + 28 x11*x19 + 2112 x12^2 + 2112 x12*x13 + 320 x12*x14 + 160 x12*x15
      + 496 x12*x16 + 248 x12*x17 + 352 x12*x18 + 176 x12*x19 + 528 x13^2
      + 160 x13*x14 + 80 x13*x15 + 248 x13*x16 + 124 x13*x17 + 176 x13*x18
      + 88 x13*x19 + 432 x14^2 + 432 x14*x15 + 512 x14*x16 + 256 x14*x17
      + 1024 x14*x18 + 512 x14*x19 + 108 x15^2 + 256 x15*x16 + 128 x15*x17
      + 512 x15*x18 + 256 x15*x19 + 3232 x16^2 + 3232 x16*x17 + 3968 x16*x18
      + 1984 x16*x19 + 808 x17^2 + 1984 x17*x18 + 992 x17*x19 + 2448 x18^2
      + 2448 x18*x19 + 612 x19^2 ]/2
Subject To
 c0: 2 x0 + x1 + 2 x2 + x3 + 2 x4 + x5 + 2 x6 + x7 + 2 x8 + x9 + 2 x10 + x11
     + 2 x12 + x13 + 2 x14 + x15 + 2 x16 + x17 + 2 x18 + x19 = 3

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

Binaries
 x0 x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18 x19
End
