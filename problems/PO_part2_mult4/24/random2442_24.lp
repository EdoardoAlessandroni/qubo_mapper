\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: 216 x0 + 108 x1 - 90 x2 - 45 x3 - 1290 x4 - 645 x5 + 252 x6 + 126 x7
      - 236 x8 - 118 x9 - 224 x10 - 112 x11 - 402 x12 - 201 x13 + 364 x14
      + 182 x15 + 226 x16 + 113 x17 - 168 x18 - 84 x19 + 218 x20 + 109 x21
      + 116 x22 + 58 x23 + [ 1144 x0^2 + 1144 x0*x1 + 832 x0*x2 + 416 x0*x3
      + 624 x0*x4 + 312 x0*x5 + 528 x0*x6 + 264 x0*x7 + 736 x0*x8 + 368 x0*x9
      + 192 x0*x10 + 96 x0*x11 + 608 x0*x12 + 304 x0*x13 + 64 x0*x14 + 32 x0*x15
      + 784 x0*x16 + 392 x0*x17 + 928 x0*x18 + 464 x0*x19 + 480 x0*x20
      + 240 x0*x21 - 192 x0*x22 - 96 x0*x23 + 286 x1^2 + 416 x1*x2 + 208 x1*x3
      + 312 x1*x4 + 156 x1*x5 + 264 x1*x6 + 132 x1*x7 + 368 x1*x8 + 184 x1*x9
      + 96 x1*x10 + 48 x1*x11 + 304 x1*x12 + 152 x1*x13 + 32 x1*x14 + 16 x1*x15
      + 392 x1*x16 + 196 x1*x17 + 464 x1*x18 + 232 x1*x19 + 240 x1*x20
      + 120 x1*x21 - 96 x1*x22 - 48 x1*x23 + 664 x2^2 + 664 x2*x3 + 96 x2*x4
      + 48 x2*x5 + 496 x2*x6 + 248 x2*x7 + 1376 x2*x8 + 688 x2*x9 + 512 x2*x10
      + 256 x2*x11 + 592 x2*x12 + 296 x2*x13 + 256 x2*x14 + 128 x2*x15
      + 1008 x2*x16 + 504 x2*x17 + 656 x2*x18 + 328 x2*x19 + 496 x2*x20
      + 248 x2*x21 + 16 x2*x22 + 8 x2*x23 + 166 x3^2 + 48 x3*x4 + 24 x3*x5
      + 248 x3*x6 + 124 x3*x7 + 688 x3*x8 + 344 x3*x9 + 256 x3*x10 + 128 x3*x11
      + 296 x3*x12 + 148 x3*x13 + 128 x3*x14 + 64 x3*x15 + 504 x3*x16
      + 252 x3*x17 + 328 x3*x18 + 164 x3*x19 + 248 x3*x20 + 124 x3*x21
      + 8 x3*x22 + 4 x3*x23 + 2080 x4^2 + 2080 x4*x5 + 816 x4*x6 + 408 x4*x7
      + 288 x4*x8 + 144 x4*x9 + 624 x4*x10 + 312 x4*x11 + 528 x4*x12
      + 264 x4*x13 + 640 x4*x14 + 320 x4*x15 + 400 x4*x18 + 200 x4*x19
      - 368 x4*x20 - 184 x4*x21 + 352 x4*x22 + 176 x4*x23 + 520 x5^2 + 408 x5*x6
      + 204 x5*x7 + 144 x5*x8 + 72 x5*x9 + 312 x5*x10 + 156 x5*x11 + 264 x5*x12
      + 132 x5*x13 + 320 x5*x14 + 160 x5*x15 + 200 x5*x18 + 100 x5*x19
      - 184 x5*x20 - 92 x5*x21 + 176 x5*x22 + 88 x5*x23 + 536 x6^2 + 536 x6*x7
      + 240 x6*x8 + 120 x6*x9 + 128 x6*x10 + 64 x6*x11 + 304 x6*x12 + 152 x6*x13
      + 304 x6*x14 + 152 x6*x15 + 256 x6*x16 + 128 x6*x17 + 224 x6*x18
      + 112 x6*x19 + 96 x6*x20 + 48 x6*x21 + 544 x6*x22 + 272 x6*x23 + 134 x7^2
      + 120 x7*x8 + 60 x7*x9 + 64 x7*x10 + 32 x7*x11 + 152 x7*x12 + 76 x7*x13
      + 152 x7*x14 + 76 x7*x15 + 128 x7*x16 + 64 x7*x17 + 112 x7*x18 + 56 x7*x19
      + 48 x7*x20 + 24 x7*x21 + 272 x7*x22 + 136 x7*x23 + 2448 x8^2 + 2448 x8*x9
      + 832 x8*x10 + 416 x8*x11 + 336 x8*x12 + 168 x8*x13 - 144 x8*x14
      - 72 x8*x15 + 1808 x8*x16 + 904 x8*x17 + 848 x8*x18 + 424 x8*x19
      + 768 x8*x20 + 384 x8*x21 + 320 x8*x22 + 160 x8*x23 + 612 x9^2
      + 416 x9*x10 + 208 x9*x11 + 168 x9*x12 + 84 x9*x13 - 72 x9*x14 - 36 x9*x15
      + 904 x9*x16 + 452 x9*x17 + 424 x9*x18 + 212 x9*x19 + 384 x9*x20
      + 192 x9*x21 + 160 x9*x22 + 80 x9*x23 + 1152 x10^2 + 1152 x10*x11
      + 352 x10*x12 + 176 x10*x13 + 464 x10*x14 + 232 x10*x15 + 480 x10*x16
      + 240 x10*x17 + 912 x10*x18 + 456 x10*x19 + 176 x10*x20 + 88 x10*x21
      + 288 x10*x22 + 144 x10*x23 + 288 x11^2 + 176 x11*x12 + 88 x11*x13
      + 232 x11*x14 + 116 x11*x15 + 240 x11*x16 + 120 x11*x17 + 456 x11*x18
      + 228 x11*x19 + 88 x11*x20 + 44 x11*x21 + 144 x11*x22 + 72 x11*x23
      + 424 x12^2 + 424 x12*x13 + 432 x12*x14 + 216 x12*x15 + 304 x12*x16
      + 152 x12*x17 + 432 x12*x18 + 216 x12*x19 + 304 x12*x20 + 152 x12*x21
      + 106 x13^2 + 216 x13*x14 + 108 x13*x15 + 152 x13*x16 + 76 x13*x17
      + 216 x13*x18 + 108 x13*x19 + 152 x13*x20 + 76 x13*x21 + 608 x14^2
      + 608 x14*x15 + 432 x14*x16 + 216 x14*x17 + 336 x14*x18 + 168 x14*x19
      - 32 x14*x20 - 16 x14*x21 + 288 x14*x22 + 144 x14*x23 + 152 x15^2
      + 216 x15*x16 + 108 x15*x17 + 168 x15*x18 + 84 x15*x19 - 16 x15*x20
      - 8 x15*x21 + 144 x15*x22 + 72 x15*x23 + 832 x16^2 + 832 x16*x17
      + 784 x16*x18 + 392 x16*x19 + 512 x16*x20 + 256 x16*x21 - 112 x16*x22
      - 56 x16*x23 + 208 x17^2 + 392 x17*x18 + 196 x17*x19 + 256 x17*x20
      + 128 x17*x21 - 56 x17*x22 - 28 x17*x23 + 456 x18^2 + 456 x18*x19
      + 352 x18*x20 + 176 x18*x21 + 16 x18*x22 + 8 x18*x23 + 114 x19^2
      + 176 x19*x20 + 88 x19*x21 + 8 x19*x22 + 4 x19*x23 + 496 x20^2
      + 496 x20*x21 - 32 x20*x22 - 16 x20*x23 + 124 x21^2 - 16 x21*x22
      - 8 x21*x23 + 560 x22^2 + 560 x22*x23 + 140 x23^2 ]/2
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
