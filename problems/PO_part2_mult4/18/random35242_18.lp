\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: - 78 x0 - 39 x1 + 308 x2 + 154 x3 + 56 x4 + 28 x5 + 562 x6 + 281 x7
      - 296 x8 - 148 x9 - 298 x10 - 149 x11 - 68 x12 - 34 x13 - 434 x14
      - 217 x15 + 400 x16 + 200 x17 + [ 976 x0^2 + 976 x0*x1 + 1136 x0*x2
      + 568 x0*x3 + 400 x0*x4 + 200 x0*x5 + 640 x0*x6 + 320 x0*x7 + 800 x0*x8
      + 400 x0*x9 + 592 x0*x10 + 296 x0*x11 + 1008 x0*x12 + 504 x0*x13
      + 976 x0*x14 + 488 x0*x15 + 1120 x0*x16 + 560 x0*x17 + 244 x1^2
      + 568 x1*x2 + 284 x1*x3 + 200 x1*x4 + 100 x1*x5 + 320 x1*x6 + 160 x1*x7
      + 400 x1*x8 + 200 x1*x9 + 296 x1*x10 + 148 x1*x11 + 504 x1*x12
      + 252 x1*x13 + 488 x1*x14 + 244 x1*x15 + 560 x1*x16 + 280 x1*x17
      + 1128 x2^2 + 1128 x2*x3 + 560 x2*x4 + 280 x2*x5 + 656 x2*x6 + 328 x2*x7
      + 912 x2*x8 + 456 x2*x9 + 912 x2*x10 + 456 x2*x11 + 752 x2*x12
      + 376 x2*x13 + 1024 x2*x14 + 512 x2*x15 + 1440 x2*x16 + 720 x2*x17
      + 282 x3^2 + 280 x3*x4 + 140 x3*x5 + 328 x3*x6 + 164 x3*x7 + 456 x3*x8
      + 228 x3*x9 + 456 x3*x10 + 228 x3*x11 + 376 x3*x12 + 188 x3*x13
      + 512 x3*x14 + 256 x3*x15 + 720 x3*x16 + 360 x3*x17 + 488 x4^2 + 488 x4*x5
      + 224 x4*x6 + 112 x4*x7 + 416 x4*x8 + 208 x4*x9 + 864 x4*x10 + 432 x4*x11
      + 448 x4*x12 + 224 x4*x13 + 992 x4*x14 + 496 x4*x15 + 1280 x4*x16
      + 640 x4*x17 + 122 x5^2 + 112 x5*x6 + 56 x5*x7 + 208 x5*x8 + 104 x5*x9
      + 432 x5*x10 + 216 x5*x11 + 224 x5*x12 + 112 x5*x13 + 496 x5*x14
      + 248 x5*x15 + 640 x5*x16 + 320 x5*x17 + 608 x6^2 + 608 x6*x7 + 624 x6*x8
      + 312 x6*x9 + 304 x6*x10 + 152 x6*x11 + 640 x6*x12 + 320 x6*x13
      + 752 x6*x14 + 376 x6*x15 + 320 x6*x16 + 160 x6*x17 + 152 x7^2 + 312 x7*x8
      + 156 x7*x9 + 152 x7*x10 + 76 x7*x11 + 320 x7*x12 + 160 x7*x13
      + 376 x7*x14 + 188 x7*x15 + 160 x7*x16 + 80 x7*x17 + 536 x8^2 + 536 x8*x9
      + 352 x8*x10 + 176 x8*x11 + 1328 x8*x12 + 664 x8*x13 + 1440 x8*x14
      + 720 x8*x15 + 320 x8*x16 + 160 x8*x17 + 134 x9^2 + 176 x9*x10 + 88 x9*x11
      + 664 x9*x12 + 332 x9*x13 + 720 x9*x14 + 360 x9*x15 + 160 x9*x16
      + 80 x9*x17 + 768 x10^2 + 768 x10*x11 + 368 x10*x12 + 184 x10*x13
      + 720 x10*x14 + 360 x10*x15 + 2368 x10*x16 + 1184 x10*x17 + 192 x11^2
      + 184 x11*x12 + 92 x11*x13 + 360 x11*x14 + 180 x11*x15 + 1184 x11*x16
      + 592 x11*x17 + 2160 x12^2 + 2160 x12*x13 + 1152 x12*x14 + 576 x12*x15
      + 224 x12*x16 + 112 x12*x17 + 540 x13^2 + 576 x13*x14 + 288 x13*x15
      + 112 x13*x16 + 56 x13*x17 + 3328 x14^2 + 3328 x14*x15 + 768 x14*x16
      + 384 x14*x17 + 832 x15^2 + 384 x15*x16 + 192 x15*x17 + 2528 x16^2
      + 2528 x16*x17 + 632 x17^2 ]/2
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
