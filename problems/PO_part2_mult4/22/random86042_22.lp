\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: - 10 x0 - 5 x1 - 262 x2 - 131 x3 - 96 x4 - 48 x5 + 40 x6 + 20 x7 + 226 x8
      + 113 x9 + 252 x10 + 126 x11 - 458 x12 - 229 x13 - 590 x14 - 295 x15
      + 278 x16 + 139 x17 - 258 x18 - 129 x19 + 90 x20 + 45 x21 + [ 912 x0^2
      + 912 x0*x1 + 320 x0*x2 + 160 x0*x3 + 1184 x0*x4 + 592 x0*x5 + 640 x0*x6
      + 320 x0*x7 + 32 x0*x8 + 16 x0*x9 + 1072 x0*x10 + 536 x0*x11 + 112 x0*x12
      + 56 x0*x13 + 448 x0*x14 + 224 x0*x15 + 80 x0*x16 + 40 x0*x17 + 256 x0*x18
      + 128 x0*x19 + 192 x0*x20 + 96 x0*x21 + 228 x1^2 + 160 x1*x2 + 80 x1*x3
      + 592 x1*x4 + 296 x1*x5 + 320 x1*x6 + 160 x1*x7 + 16 x1*x8 + 8 x1*x9
      + 536 x1*x10 + 268 x1*x11 + 56 x1*x12 + 28 x1*x13 + 224 x1*x14
      + 112 x1*x15 + 40 x1*x16 + 20 x1*x17 + 128 x1*x18 + 64 x1*x19 + 96 x1*x20
      + 48 x1*x21 + 696 x2^2 + 696 x2*x3 + 400 x2*x4 + 200 x2*x5 + 1280 x2*x6
      + 640 x2*x7 + 960 x2*x8 + 480 x2*x9 + 384 x2*x10 + 192 x2*x11 + 544 x2*x12
      + 272 x2*x13 + 672 x2*x14 + 336 x2*x15 + 688 x2*x16 + 344 x2*x17
      + 368 x2*x18 + 184 x2*x19 + 1248 x2*x20 + 624 x2*x21 + 174 x3^2
      + 200 x3*x4 + 100 x3*x5 + 640 x3*x6 + 320 x3*x7 + 480 x3*x8 + 240 x3*x9
      + 192 x3*x10 + 96 x3*x11 + 272 x3*x12 + 136 x3*x13 + 336 x3*x14
      + 168 x3*x15 + 344 x3*x16 + 172 x3*x17 + 184 x3*x18 + 92 x3*x19
      + 624 x3*x20 + 312 x3*x21 + 512 x4^2 + 512 x4*x5 + 736 x4*x6 + 368 x4*x7
      + 144 x4*x8 + 72 x4*x9 + 592 x4*x10 + 296 x4*x11 + 64 x4*x12 + 32 x4*x13
      + 368 x4*x14 + 184 x4*x15 + 320 x4*x16 + 160 x4*x17 + 96 x4*x18
      + 48 x4*x19 + 576 x4*x20 + 288 x4*x21 + 128 x5^2 + 368 x5*x6 + 184 x5*x7
      + 72 x5*x8 + 36 x5*x9 + 296 x5*x10 + 148 x5*x11 + 32 x5*x12 + 16 x5*x13
      + 184 x5*x14 + 92 x5*x15 + 160 x5*x16 + 80 x5*x17 + 48 x5*x18 + 24 x5*x19
      + 288 x5*x20 + 144 x5*x21 + 984 x6^2 + 984 x6*x7 + 1232 x6*x8 + 616 x6*x9
      + 496 x6*x10 + 248 x6*x11 + 416 x6*x12 + 208 x6*x13 + 768 x6*x14
      + 384 x6*x15 + 784 x6*x16 + 392 x6*x17 - 48 x6*x18 - 24 x6*x19
      + 1728 x6*x20 + 864 x6*x21 + 246 x7^2 + 616 x7*x8 + 308 x7*x9 + 248 x7*x10
      + 124 x7*x11 + 208 x7*x12 + 104 x7*x13 + 384 x7*x14 + 192 x7*x15
      + 392 x7*x16 + 196 x7*x17 - 24 x7*x18 - 12 x7*x19 + 864 x7*x20
      + 432 x7*x21 + 832 x8^2 + 832 x8*x9 + 256 x8*x10 + 128 x8*x11 + 80 x8*x12
      + 40 x8*x13 + 272 x8*x14 + 136 x8*x15 + 416 x8*x16 + 208 x8*x17
      - 416 x8*x18 - 208 x8*x19 + 1712 x8*x20 + 856 x8*x21 + 208 x9^2
      + 128 x9*x10 + 64 x9*x11 + 40 x9*x12 + 20 x9*x13 + 136 x9*x14 + 68 x9*x15
      + 208 x9*x16 + 104 x9*x17 - 208 x9*x18 - 104 x9*x19 + 856 x9*x20
      + 428 x9*x21 + 536 x10^2 + 536 x10*x11 - 80 x10*x12 - 40 x10*x13
      + 288 x10*x14 + 144 x10*x15 + 80 x10*x16 + 40 x10*x17 + 224 x10*x18
      + 112 x10*x19 - 256 x10*x20 - 128 x10*x21 + 134 x11^2 - 40 x11*x12
      - 20 x11*x13 + 144 x11*x14 + 72 x11*x15 + 40 x11*x16 + 20 x11*x17
      + 112 x11*x18 + 56 x11*x19 - 128 x11*x20 - 64 x11*x21 + 816 x12^2
      + 816 x12*x13 + 592 x12*x14 + 296 x12*x15 + 272 x12*x16 + 136 x12*x17
      + 576 x12*x18 + 288 x12*x19 + 288 x12*x20 + 144 x12*x21 + 204 x13^2
      + 296 x13*x14 + 148 x13*x15 + 136 x13*x16 + 68 x13*x17 + 288 x13*x18
      + 144 x13*x19 + 144 x13*x20 + 72 x13*x21 + 392 x14^2 + 392 x14*x15
      + 432 x14*x16 + 216 x14*x17 + 320 x14*x18 + 160 x14*x19 + 368 x14*x20
      + 184 x14*x21 + 98 x15^2 + 216 x15*x16 + 108 x15*x17 + 160 x15*x18
      + 80 x15*x19 + 184 x15*x20 + 92 x15*x21 + 600 x16^2 + 600 x16*x17
      + 208 x16*x18 + 104 x16*x19 + 1088 x16*x20 + 544 x16*x21 + 150 x17^2
      + 104 x17*x18 + 52 x17*x19 + 544 x17*x20 + 272 x17*x21 + 808 x18^2
      + 808 x18*x19 + 480 x18*x20 + 240 x18*x21 + 202 x19^2 + 240 x19*x20
      + 120 x19*x21 + 2528 x20^2 + 2528 x20*x21 + 632 x21^2 ]/2
Subject To
 c0: 2 x0 + x1 + 2 x2 + x3 + 2 x4 + x5 + 2 x6 + x7 + 2 x8 + x9 + 2 x10 + x11
     + 2 x12 + x13 + 2 x14 + x15 + 2 x16 + x17 + 2 x18 + x19 + 2 x20 + x21 = 3

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

Binaries
 x0 x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18 x19 x20 x21
End
