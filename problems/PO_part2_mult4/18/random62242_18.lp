\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: 284 x0 + 142 x1 - 90 x2 - 45 x3 + 42 x4 + 21 x5 - 616 x6 - 308 x7 - 254 x8
      - 127 x9 + 178 x10 + 89 x11 + 328 x12 + 164 x13 - 1290 x14 - 645 x15
      - 420 x16 - 210 x17 + [ 1104 x0^2 + 1104 x0*x1 + 1216 x0*x2 + 608 x0*x3
      + 608 x0*x4 + 304 x0*x5 + 1152 x0*x6 + 576 x0*x7 + 592 x0*x8 + 296 x0*x9
      + 1808 x0*x10 + 904 x0*x11 + 464 x0*x12 + 232 x0*x13 + 224 x0*x14
      + 112 x0*x15 + 640 x0*x16 + 320 x0*x17 + 276 x1^2 + 608 x1*x2 + 304 x1*x3
      + 304 x1*x4 + 152 x1*x5 + 576 x1*x6 + 288 x1*x7 + 296 x1*x8 + 148 x1*x9
      + 904 x1*x10 + 452 x1*x11 + 232 x1*x12 + 116 x1*x13 + 112 x1*x14
      + 56 x1*x15 + 320 x1*x16 + 160 x1*x17 + 664 x2^2 + 664 x2*x3 + 768 x2*x4
      + 384 x2*x5 + 848 x2*x6 + 424 x2*x7 + 416 x2*x8 + 208 x2*x9 + 960 x2*x10
      + 480 x2*x11 + 336 x2*x12 + 168 x2*x13 + 96 x2*x14 + 48 x2*x15
      + 320 x2*x16 + 160 x2*x17 + 166 x3^2 + 384 x3*x4 + 192 x3*x5 + 424 x3*x6
      + 212 x3*x7 + 208 x3*x8 + 104 x3*x9 + 480 x3*x10 + 240 x3*x11 + 168 x3*x12
      + 84 x3*x13 + 48 x3*x14 + 24 x3*x15 + 160 x3*x16 + 80 x3*x17 + 1144 x4^2
      + 1144 x4*x5 + 448 x4*x6 + 224 x4*x7 + 208 x4*x8 + 104 x4*x9 + 448 x4*x10
      + 224 x4*x11 - 48 x4*x12 - 24 x4*x13 - 64 x4*x14 - 32 x4*x15 + 144 x4*x16
      + 72 x4*x17 + 286 x5^2 + 224 x5*x6 + 112 x5*x7 + 104 x5*x8 + 52 x5*x9
      + 224 x5*x10 + 112 x5*x11 - 24 x5*x12 - 12 x5*x13 - 32 x5*x14 - 16 x5*x15
      + 72 x5*x16 + 36 x5*x17 + 992 x6^2 + 992 x6*x7 + 608 x6*x8 + 304 x6*x9
      + 1104 x6*x10 + 552 x6*x11 + 432 x6*x12 + 216 x6*x13 - 272 x6*x14
      - 136 x6*x15 + 512 x6*x16 + 256 x6*x17 + 248 x7^2 + 304 x7*x8 + 152 x7*x9
      + 552 x7*x10 + 276 x7*x11 + 216 x7*x12 + 108 x7*x13 - 136 x7*x14
      - 68 x7*x15 + 256 x7*x16 + 128 x7*x17 + 176 x8^2 + 176 x8*x9 + 480 x8*x10
      + 240 x8*x11 + 192 x8*x12 + 96 x8*x13 + 32 x8*x14 + 16 x8*x15 + 240 x8*x16
      + 120 x8*x17 + 44 x9^2 + 240 x9*x10 + 120 x9*x11 + 96 x9*x12 + 48 x9*x13
      + 16 x9*x14 + 8 x9*x15 + 120 x9*x16 + 60 x9*x17 + 1064 x10^2
      + 1064 x10*x11 + 384 x10*x12 + 192 x10*x13 + 400 x10*x16 + 200 x10*x17
      + 266 x11^2 + 192 x11*x12 + 96 x11*x13 + 200 x11*x16 + 100 x11*x17
      + 1128 x12^2 + 1128 x12*x13 + 864 x12*x14 + 432 x12*x15 + 400 x12*x16
      + 200 x12*x17 + 282 x13^2 + 432 x13*x14 + 216 x13*x15 + 200 x13*x16
      + 100 x13*x17 + 2080 x14^2 + 2080 x14*x15 + 544 x14*x16 + 272 x14*x17
      + 520 x15^2 + 272 x15*x16 + 136 x15*x17 + 376 x16^2 + 376 x16*x17
      + 94 x17^2 ]/2
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
