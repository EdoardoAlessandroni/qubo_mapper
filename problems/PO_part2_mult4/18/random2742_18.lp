\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: 42 x0 + 21 x1 + 376 x2 + 188 x3 - 66 x4 - 33 x5 - 288 x6 - 144 x7 - 312 x8
      - 156 x9 + 90 x10 + 45 x11 + 18 x12 + 9 x13 - 68 x14 - 34 x15 - 76 x16
      - 38 x17 + [ 1016 x0^2 + 1016 x0*x1 + 688 x0*x2 + 344 x0*x3 + 272 x0*x4
      + 136 x0*x5 - 240 x0*x6 - 120 x0*x7 + 16 x0*x8 + 8 x0*x9 + 704 x0*x10
      + 352 x0*x11 + 1008 x0*x12 + 504 x0*x13 + 1824 x0*x14 + 912 x0*x15
      + 896 x0*x16 + 448 x0*x17 + 254 x1^2 + 344 x1*x2 + 172 x1*x3 + 136 x1*x4
      + 68 x1*x5 - 120 x1*x6 - 60 x1*x7 + 8 x1*x8 + 4 x1*x9 + 352 x1*x10
      + 176 x1*x11 + 504 x1*x12 + 252 x1*x13 + 912 x1*x14 + 456 x1*x15
      + 448 x1*x16 + 224 x1*x17 + 792 x2^2 + 792 x2*x3 + 240 x2*x4 + 120 x2*x5
      - 176 x2*x6 - 88 x2*x7 + 736 x2*x10 + 368 x2*x11 + 400 x2*x12 + 200 x2*x13
      + 1584 x2*x14 + 792 x2*x15 + 336 x2*x16 + 168 x2*x17 + 198 x3^2
      + 120 x3*x4 + 60 x3*x5 - 88 x3*x6 - 44 x3*x7 + 368 x3*x10 + 184 x3*x11
      + 200 x3*x12 + 100 x3*x13 + 792 x3*x14 + 396 x3*x15 + 168 x3*x16
      + 84 x3*x17 + 296 x4^2 + 296 x4*x5 - 160 x4*x6 - 80 x4*x7 - 304 x4*x8
      - 152 x4*x9 + 176 x4*x10 + 88 x4*x11 + 464 x4*x12 + 232 x4*x13
      + 144 x4*x14 + 72 x4*x15 + 496 x4*x16 + 248 x4*x17 + 74 x5^2 - 80 x5*x6
      - 40 x5*x7 - 152 x5*x8 - 76 x5*x9 + 88 x5*x10 + 44 x5*x11 + 232 x5*x12
      + 116 x5*x13 + 72 x5*x14 + 36 x5*x15 + 248 x5*x16 + 124 x5*x17 + 320 x6^2
      + 320 x6*x7 + 80 x6*x8 + 40 x6*x9 + 80 x6*x10 + 40 x6*x11 + 208 x6*x12
      + 104 x6*x13 - 608 x6*x14 - 304 x6*x15 - 144 x6*x16 - 72 x6*x17 + 80 x7^2
      + 40 x7*x8 + 20 x7*x9 + 40 x7*x10 + 20 x7*x11 + 104 x7*x12 + 52 x7*x13
      - 304 x7*x14 - 152 x7*x15 - 72 x7*x16 - 36 x7*x17 + 1168 x8^2 + 1168 x8*x9
      + 32 x8*x10 + 16 x8*x11 - 432 x8*x12 - 216 x8*x13 - 480 x8*x14
      - 240 x8*x15 - 192 x8*x16 - 96 x8*x17 + 292 x9^2 + 16 x9*x10 + 8 x9*x11
      - 216 x9*x12 - 108 x9*x13 - 240 x9*x14 - 120 x9*x15 - 96 x9*x16
      - 48 x9*x17 + 664 x10^2 + 664 x10*x11 + 576 x10*x12 + 288 x10*x13
      + 832 x10*x14 + 416 x10*x15 + 160 x10*x16 + 80 x10*x17 + 166 x11^2
      + 288 x11*x12 + 144 x11*x13 + 416 x11*x14 + 208 x11*x15 + 80 x11*x16
      + 40 x11*x17 + 832 x12^2 + 832 x12*x13 + 384 x12*x14 + 192 x12*x15
      + 816 x12*x16 + 408 x12*x17 + 208 x13^2 + 192 x13*x14 + 96 x13*x15
      + 408 x13*x16 + 204 x13*x17 + 2160 x14^2 + 2160 x14*x15 + 192 x14*x16
      + 96 x14*x17 + 540 x15^2 + 96 x15*x16 + 48 x15*x17 + 560 x16^2
      + 560 x16*x17 + 140 x17^2 ]/2
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
