\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: 686 x0 + 343 x1 - 422 x2 - 211 x3 - 672 x4 - 336 x5 - 246 x6 - 123 x7
      - 468 x8 - 234 x9 + 90 x10 + 45 x11 + 178 x12 + 89 x13 - 1290 x14
      - 645 x15 - 114 x16 - 57 x17 + [ 1352 x0^2 + 1352 x0*x1 + 1056 x0*x2
      + 528 x0*x3 + 800 x0*x4 + 400 x0*x5 + 880 x0*x6 + 440 x0*x7 + 272 x0*x8
      + 136 x0*x9 + 1696 x0*x10 + 848 x0*x11 + 1104 x0*x12 + 552 x0*x13
      + 1120 x0*x14 + 560 x0*x15 + 480 x0*x16 + 240 x0*x17 + 338 x1^2
      + 528 x1*x2 + 264 x1*x3 + 400 x1*x4 + 200 x1*x5 + 440 x1*x6 + 220 x1*x7
      + 136 x1*x8 + 68 x1*x9 + 848 x1*x10 + 424 x1*x11 + 552 x1*x12 + 276 x1*x13
      + 560 x1*x14 + 280 x1*x15 + 240 x1*x16 + 120 x1*x17 + 576 x2^2 + 576 x2*x3
      + 448 x2*x4 + 224 x2*x5 + 592 x2*x6 + 296 x2*x7 + 624 x2*x8 + 312 x2*x9
      + 1168 x2*x10 + 584 x2*x11 + 672 x2*x12 + 336 x2*x13 + 432 x2*x14
      + 216 x2*x15 + 304 x2*x16 + 152 x2*x17 + 144 x3^2 + 224 x3*x4 + 112 x3*x5
      + 296 x3*x6 + 148 x3*x7 + 312 x3*x8 + 156 x3*x9 + 584 x3*x10 + 292 x3*x11
      + 336 x3*x12 + 168 x3*x13 + 216 x3*x14 + 108 x3*x15 + 152 x3*x16
      + 76 x3*x17 + 1048 x4^2 + 1048 x4*x5 + 1024 x4*x6 + 512 x4*x7 + 128 x4*x8
      + 64 x4*x9 + 1600 x4*x10 + 800 x4*x11 + 480 x4*x12 + 240 x4*x13
      + 640 x4*x14 + 320 x4*x15 - 32 x4*x16 - 16 x4*x17 + 262 x5^2 + 512 x5*x6
      + 256 x5*x7 + 64 x5*x8 + 32 x5*x9 + 800 x5*x10 + 400 x5*x11 + 240 x5*x12
      + 120 x5*x13 + 320 x5*x14 + 160 x5*x15 - 16 x5*x16 - 8 x5*x17 + 688 x6^2
      + 688 x6*x7 + 64 x6*x8 + 32 x6*x9 + 1808 x6*x10 + 904 x6*x11 + 1040 x6*x12
      + 520 x6*x13 + 480 x6*x14 + 240 x6*x15 - 16 x6*x16 - 8 x6*x17 + 172 x7^2
      + 32 x7*x8 + 16 x7*x9 + 904 x7*x10 + 452 x7*x11 + 520 x7*x12 + 260 x7*x13
      + 240 x7*x14 + 120 x7*x15 - 8 x7*x16 - 4 x7*x17 + 504 x8^2 + 504 x8*x9
      - 240 x8*x10 - 120 x8*x11 + 432 x8*x12 + 216 x8*x13 - 192 x8*x14
      - 96 x8*x15 + 32 x8*x16 + 16 x8*x17 + 126 x9^2 - 120 x9*x10 - 60 x9*x11
      + 216 x9*x12 + 108 x9*x13 - 96 x9*x14 - 48 x9*x15 + 16 x9*x16 + 8 x9*x17
      + 2528 x10^2 + 2528 x10*x11 + 1632 x10*x12 + 816 x10*x13 + 128 x10*x14
      + 64 x10*x15 + 368 x10*x16 + 184 x10*x17 + 632 x11^2 + 816 x11*x12
      + 408 x11*x13 + 64 x11*x14 + 32 x11*x15 + 184 x11*x16 + 92 x11*x17
      + 1064 x12^2 + 1064 x12*x13 + 128 x12*x16 + 64 x12*x17 + 266 x13^2
      + 64 x13*x16 + 32 x13*x17 + 2080 x14^2 + 2080 x14*x15 + 352 x14*x16
      + 176 x14*x17 + 520 x15^2 + 176 x15*x16 + 88 x15*x17 + 328 x16^2
      + 328 x16*x17 + 82 x17^2 ]/2
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
