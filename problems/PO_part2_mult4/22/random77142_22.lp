\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: 42 x0 + 21 x1 + 586 x2 + 293 x3 - 52 x4 - 26 x5 - 78 x6 - 39 x7 - 290 x8
      - 145 x9 - 110 x10 - 55 x11 + 584 x12 + 292 x13 - 242 x14 - 121 x15
      - 590 x16 - 295 x17 + 116 x18 + 58 x19 - 44 x20 - 22 x21 + [ 1144 x0^2
      + 1144 x0*x1 - 480 x0*x2 - 240 x0*x3 + 224 x0*x4 + 112 x0*x5 + 992 x0*x6
      + 496 x0*x7 + 576 x0*x8 + 288 x0*x9 + 1056 x0*x10 + 528 x0*x11
      + 416 x0*x12 + 208 x0*x13 + 1040 x0*x14 + 520 x0*x15 + 656 x0*x16
      + 328 x0*x17 + 384 x0*x18 + 192 x0*x19 + 752 x0*x20 + 376 x0*x21
      + 286 x1^2 - 240 x1*x2 - 120 x1*x3 + 112 x1*x4 + 56 x1*x5 + 496 x1*x6
      + 248 x1*x7 + 288 x1*x8 + 144 x1*x9 + 528 x1*x10 + 264 x1*x11 + 208 x1*x12
      + 104 x1*x13 + 520 x1*x14 + 260 x1*x15 + 328 x1*x16 + 164 x1*x17
      + 192 x1*x18 + 96 x1*x19 + 376 x1*x20 + 188 x1*x21 + 1240 x2^2
      + 1240 x2*x3 + 336 x2*x4 + 168 x2*x5 + 224 x2*x6 + 112 x2*x7 + 16 x2*x8
      + 8 x2*x9 + 720 x2*x10 + 360 x2*x11 - 80 x2*x12 - 40 x2*x13 + 16 x2*x14
      + 8 x2*x15 + 384 x2*x16 + 192 x2*x17 - 160 x2*x18 - 80 x2*x19 + 656 x2*x20
      + 328 x2*x21 + 310 x3^2 + 168 x3*x4 + 84 x3*x5 + 112 x3*x6 + 56 x3*x7
      + 8 x3*x8 + 4 x3*x9 + 360 x3*x10 + 180 x3*x11 - 40 x3*x12 - 20 x3*x13
      + 8 x3*x14 + 4 x3*x15 + 192 x3*x16 + 96 x3*x17 - 80 x3*x18 - 40 x3*x19
      + 328 x3*x20 + 164 x3*x21 + 376 x4^2 + 376 x4*x5 + 384 x4*x6 + 192 x4*x7
      + 528 x4*x10 + 264 x4*x11 + 272 x4*x12 + 136 x4*x13 + 352 x4*x14
      + 176 x4*x15 + 208 x4*x16 + 104 x4*x17 - 208 x4*x18 - 104 x4*x19
      - 224 x4*x20 - 112 x4*x21 + 94 x5^2 + 192 x5*x6 + 96 x5*x7 + 264 x5*x10
      + 132 x5*x11 + 136 x5*x12 + 68 x5*x13 + 176 x5*x14 + 88 x5*x15
      + 104 x5*x16 + 52 x5*x17 - 104 x5*x18 - 52 x5*x19 - 112 x5*x20 - 56 x5*x21
      + 976 x6^2 + 976 x6*x7 + 576 x6*x8 + 288 x6*x9 + 912 x6*x10 + 456 x6*x11
      + 224 x6*x12 + 112 x6*x13 + 912 x6*x14 + 456 x6*x15 + 608 x6*x16
      + 304 x6*x17 - 160 x6*x18 - 80 x6*x19 + 512 x6*x20 + 256 x6*x21 + 244 x7^2
      + 288 x7*x8 + 144 x7*x9 + 456 x7*x10 + 228 x7*x11 + 112 x7*x12 + 56 x7*x13
      + 456 x7*x14 + 228 x7*x15 + 304 x7*x16 + 152 x7*x17 - 80 x7*x18
      - 40 x7*x19 + 256 x7*x20 + 128 x7*x21 + 752 x8^2 + 752 x8*x9 + 1088 x8*x10
      + 544 x8*x11 + 128 x8*x12 + 64 x8*x13 + 560 x8*x14 + 280 x8*x15
      + 544 x8*x16 + 272 x8*x17 + 176 x8*x18 + 88 x8*x19 + 1200 x8*x20
      + 600 x8*x21 + 188 x9^2 + 544 x9*x10 + 272 x9*x11 + 64 x9*x12 + 32 x9*x13
      + 280 x9*x14 + 140 x9*x15 + 272 x9*x16 + 136 x9*x17 + 88 x9*x18
      + 44 x9*x19 + 600 x9*x20 + 300 x9*x21 + 1224 x10^2 + 1224 x10*x11
      + 128 x10*x12 + 64 x10*x13 + 1008 x10*x14 + 504 x10*x15 + 1120 x10*x16
      + 560 x10*x17 + 32 x10*x18 + 16 x10*x19 + 800 x10*x20 + 400 x10*x21
      + 306 x11^2 + 64 x11*x12 + 32 x11*x13 + 504 x11*x14 + 252 x11*x15
      + 560 x11*x16 + 280 x11*x17 + 16 x11*x18 + 8 x11*x19 + 400 x11*x20
      + 200 x11*x21 + 1304 x12^2 + 1304 x12*x13 + 896 x12*x14 + 448 x12*x15
      + 96 x12*x16 + 48 x12*x17 - 32 x12*x18 - 16 x12*x19 - 192 x12*x20
      - 96 x12*x21 + 326 x13^2 + 448 x13*x14 + 224 x13*x15 + 48 x13*x16
      + 24 x13*x17 - 16 x13*x18 - 8 x13*x19 - 96 x13*x20 - 48 x13*x21
      + 672 x14^2 + 672 x14*x15 + 592 x14*x16 + 296 x14*x17 + 112 x14*x18
      + 56 x14*x19 + 400 x14*x20 + 200 x14*x21 + 168 x15^2 + 296 x15*x16
      + 148 x15*x17 + 56 x15*x18 + 28 x15*x19 + 200 x15*x20 + 100 x15*x21
      + 392 x16^2 + 392 x16*x17 - 80 x16*x18 - 40 x16*x19 + 400 x16*x20
      + 200 x16*x21 + 98 x17^2 - 40 x17*x18 - 20 x17*x19 + 200 x17*x20
      + 100 x17*x21 + 560 x18^2 + 560 x18*x19 + 288 x18*x20 + 144 x18*x21
      + 140 x19^2 + 144 x19*x20 + 72 x19*x21 + 1168 x20^2 + 1168 x20*x21
      + 292 x21^2 ]/2
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
