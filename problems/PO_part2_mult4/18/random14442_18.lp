\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: - 590 x0 - 295 x1 - 236 x2 - 118 x3 - 78 x4 - 39 x5 + 54 x6 + 27 x7
      + 90 x8 + 45 x9 + 328 x10 + 164 x11 - 340 x12 - 170 x13 + 138 x14 + 69 x15
      + 870 x16 + 435 x17 + [ 392 x0^2 + 392 x0*x1 + 144 x0*x2 + 72 x0*x3
      + 608 x0*x4 + 304 x0*x5 + 320 x0*x6 + 160 x0*x7 + 432 x0*x8 + 216 x0*x9
      - 96 x0*x10 - 48 x0*x11 + 320 x0*x12 + 160 x0*x13 + 144 x0*x14 + 72 x0*x15
      + 592 x0*x16 + 296 x0*x17 + 98 x1^2 + 72 x1*x2 + 36 x1*x3 + 304 x1*x4
      + 152 x1*x5 + 160 x1*x6 + 80 x1*x7 + 216 x1*x8 + 108 x1*x9 - 48 x1*x10
      - 24 x1*x11 + 160 x1*x12 + 80 x1*x13 + 72 x1*x14 + 36 x1*x15 + 296 x1*x16
      + 148 x1*x17 + 464 x2^2 + 464 x2*x3 + 480 x2*x4 + 240 x2*x5 + 464 x2*x6
      + 232 x2*x7 + 64 x2*x8 + 32 x2*x9 - 16 x2*x10 - 8 x2*x11 + 448 x2*x12
      + 224 x2*x13 - 96 x2*x14 - 48 x2*x15 + 96 x2*x16 + 48 x2*x17 + 116 x3^2
      + 240 x3*x4 + 120 x3*x5 + 232 x3*x6 + 116 x3*x7 + 32 x3*x8 + 16 x3*x9
      - 8 x3*x10 - 4 x3*x11 + 224 x3*x12 + 112 x3*x13 - 48 x3*x14 - 24 x3*x15
      + 48 x3*x16 + 24 x3*x17 + 976 x4^2 + 976 x4*x5 + 320 x4*x6 + 160 x4*x7
      - 96 x4*x8 - 48 x4*x9 + 480 x4*x10 + 240 x4*x11 + 640 x4*x12 + 320 x4*x13
      - 160 x4*x14 - 80 x4*x15 + 496 x4*x16 + 248 x4*x17 + 244 x5^2 + 160 x5*x6
      + 80 x5*x7 - 48 x5*x8 - 24 x5*x9 + 240 x5*x10 + 120 x5*x11 + 320 x5*x12
      + 160 x5*x13 - 80 x5*x14 - 40 x5*x15 + 248 x5*x16 + 124 x5*x17 + 392 x6^2
      + 392 x6*x7 + 192 x6*x8 + 96 x6*x9 - 496 x6*x10 - 248 x6*x11 + 288 x6*x12
      + 144 x6*x13 - 32 x6*x14 - 16 x6*x15 + 80 x6*x16 + 40 x6*x17 + 98 x7^2
      + 96 x7*x8 + 48 x7*x9 - 248 x7*x10 - 124 x7*x11 + 144 x7*x12 + 72 x7*x13
      - 16 x7*x14 - 8 x7*x15 + 40 x7*x16 + 20 x7*x17 + 664 x8^2 + 664 x8*x9
      + 736 x8*x10 + 368 x8*x11 + 368 x8*x12 + 184 x8*x13 + 160 x8*x14
      + 80 x8*x15 + 1056 x8*x16 + 528 x8*x17 + 166 x9^2 + 368 x9*x10
      + 184 x9*x11 + 184 x9*x12 + 92 x9*x13 + 80 x9*x14 + 40 x9*x15 + 528 x9*x16
      + 264 x9*x17 + 1128 x10^2 + 1128 x10*x11 + 96 x10*x12 + 48 x10*x13
      + 48 x10*x14 + 24 x10*x15 + 1152 x10*x16 + 576 x10*x17 + 282 x11^2
      + 48 x11*x12 + 24 x11*x13 + 24 x11*x14 + 12 x11*x15 + 576 x11*x16
      + 288 x11*x17 + 528 x12^2 + 528 x12*x13 + 32 x12*x14 + 16 x12*x15
      + 208 x12*x16 + 104 x12*x17 + 132 x13^2 + 16 x13*x14 + 8 x13*x15
      + 104 x13*x16 + 52 x13*x17 + 1128 x14^2 + 1128 x14*x15 + 272 x14*x16
      + 136 x14*x17 + 282 x15^2 + 136 x15*x16 + 68 x15*x17 + 1032 x16^2
      + 1032 x16*x17 + 258 x17^2 ]/2
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
