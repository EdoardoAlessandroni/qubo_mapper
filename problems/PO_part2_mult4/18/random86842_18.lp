\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: 200 x0 + 100 x1 + 870 x2 + 435 x3 - 128 x4 - 64 x5 + 400 x6 + 200 x7
      + 628 x8 + 314 x9 - 66 x10 - 33 x11 - 298 x12 - 149 x13 + 172 x14 + 86 x15
      - 88 x16 - 44 x17 + [ 808 x0^2 + 808 x0*x1 + 368 x0*x2 + 184 x0*x3
      + 80 x0*x4 + 40 x0*x5 + 1456 x0*x6 + 728 x0*x7 + 512 x0*x8 + 256 x0*x9
      + 272 x0*x10 + 136 x0*x11 + 880 x0*x12 + 440 x0*x13 + 640 x0*x14
      + 320 x0*x15 + 448 x0*x16 + 224 x0*x17 + 202 x1^2 + 184 x1*x2 + 92 x1*x3
      + 40 x1*x4 + 20 x1*x5 + 728 x1*x6 + 364 x1*x7 + 256 x1*x8 + 128 x1*x9
      + 136 x1*x10 + 68 x1*x11 + 440 x1*x12 + 220 x1*x13 + 320 x1*x14
      + 160 x1*x15 + 224 x1*x16 + 112 x1*x17 + 1032 x2^2 + 1032 x2*x3
      + 240 x2*x4 + 120 x2*x5 + 1024 x2*x6 + 512 x2*x7 + 512 x2*x8 + 256 x2*x9
      + 448 x2*x10 + 224 x2*x11 + 512 x2*x12 + 256 x2*x13 + 992 x2*x14
      + 496 x2*x15 + 448 x2*x16 + 224 x2*x17 + 258 x3^2 + 120 x3*x4 + 60 x3*x5
      + 512 x3*x6 + 256 x3*x7 + 256 x3*x8 + 128 x3*x9 + 224 x3*x10 + 112 x3*x11
      + 256 x3*x12 + 128 x3*x13 + 496 x3*x14 + 248 x3*x15 + 224 x3*x16
      + 112 x3*x17 + 368 x4^2 + 368 x4*x5 + 496 x4*x6 + 248 x4*x7 - 96 x4*x8
      - 48 x4*x9 - 16 x4*x10 - 8 x4*x11 + 336 x4*x12 + 168 x4*x13 + 48 x4*x14
      + 24 x4*x15 + 384 x4*x16 + 192 x4*x17 + 92 x5^2 + 248 x5*x6 + 124 x5*x7
      - 48 x5*x8 - 24 x5*x9 - 8 x5*x10 - 4 x5*x11 + 168 x5*x12 + 84 x5*x13
      + 24 x5*x14 + 12 x5*x15 + 192 x5*x16 + 96 x5*x17 + 2528 x6^2 + 2528 x6*x7
      + 1904 x6*x8 + 952 x6*x9 + 992 x6*x10 + 496 x6*x11 + 2368 x6*x12
      + 1184 x6*x13 - 128 x6*x14 - 64 x6*x15 + 1168 x6*x16 + 584 x6*x17
      + 632 x7^2 + 952 x7*x8 + 476 x7*x9 + 496 x7*x10 + 248 x7*x11 + 1184 x7*x12
      + 592 x7*x13 - 64 x7*x14 - 32 x7*x15 + 584 x7*x16 + 292 x7*x17 + 2112 x8^2
      + 2112 x8*x9 + 512 x8*x10 + 256 x8*x11 + 1056 x8*x12 + 528 x8*x13
      + 1168 x8*x14 + 584 x8*x15 + 208 x8*x16 + 104 x8*x17 + 528 x9^2
      + 256 x9*x10 + 128 x9*x11 + 528 x9*x12 + 264 x9*x13 + 584 x9*x14
      + 292 x9*x15 + 104 x9*x16 + 52 x9*x17 + 296 x10^2 + 296 x10*x11
      + 352 x10*x12 + 176 x10*x13 + 352 x10*x14 + 176 x10*x15 + 160 x10*x16
      + 80 x10*x17 + 74 x11^2 + 176 x11*x12 + 88 x11*x13 + 176 x11*x14
      + 88 x11*x15 + 80 x11*x16 + 40 x11*x17 + 768 x12^2 + 768 x12*x13
      + 128 x12*x14 + 64 x12*x15 + 752 x12*x16 + 376 x12*x17 + 192 x13^2
      + 64 x13*x14 + 32 x13*x15 + 376 x13*x16 + 188 x13*x17 + 1072 x14^2
      + 1072 x14*x15 + 160 x14*x16 + 80 x14*x17 + 268 x15^2 + 80 x15*x16
      + 40 x15*x17 + 472 x16^2 + 472 x16*x17 + 118 x17^2 ]/2
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
