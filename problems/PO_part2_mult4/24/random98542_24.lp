\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: - 200 x0 - 100 x1 + 226 x2 + 113 x3 - 402 x4 - 201 x5 + 584 x6 + 292 x7
      + 70 x8 + 35 x9 - 156 x10 - 78 x11 + 280 x12 + 140 x13 + 72 x14 + 36 x15
      + 36 x16 + 18 x17 + 364 x18 + 182 x19 - 258 x20 - 129 x21 + 10 x22 + 5 x23
      + [ 248 x0^2 + 248 x0*x1 - 224 x0*x2 - 112 x0*x3 + 160 x0*x4 + 80 x0*x5
      - 224 x0*x6 - 112 x0*x7 + 144 x0*x8 + 72 x0*x9 + 96 x0*x10 + 48 x0*x11
      + 16 x0*x12 + 8 x0*x13 - 112 x0*x14 - 56 x0*x15 + 192 x0*x16 + 96 x0*x17
      + 112 x0*x18 + 56 x0*x19 + 240 x0*x20 + 120 x0*x21 + 448 x0*x22
      + 224 x0*x23 + 62 x1^2 - 112 x1*x2 - 56 x1*x3 + 80 x1*x4 + 40 x1*x5
      - 112 x1*x6 - 56 x1*x7 + 72 x1*x8 + 36 x1*x9 + 48 x1*x10 + 24 x1*x11
      + 8 x1*x12 + 4 x1*x13 - 56 x1*x14 - 28 x1*x15 + 96 x1*x16 + 48 x1*x17
      + 56 x1*x18 + 28 x1*x19 + 120 x1*x20 + 60 x1*x21 + 224 x1*x22 + 112 x1*x23
      + 832 x2^2 + 832 x2*x3 + 304 x2*x4 + 152 x2*x5 + 928 x2*x6 + 464 x2*x7
      + 160 x2*x8 + 80 x2*x9 + 912 x2*x10 + 456 x2*x11 + 1040 x2*x12
      + 520 x2*x13 + 640 x2*x14 + 320 x2*x15 + 944 x2*x16 + 472 x2*x17
      + 432 x2*x18 + 216 x2*x19 - 416 x2*x20 - 208 x2*x21 - 64 x2*x22
      - 32 x2*x23 + 208 x3^2 + 152 x3*x4 + 76 x3*x5 + 464 x3*x6 + 232 x3*x7
      + 80 x3*x8 + 40 x3*x9 + 456 x3*x10 + 228 x3*x11 + 520 x3*x12 + 260 x3*x13
      + 320 x3*x14 + 160 x3*x15 + 472 x3*x16 + 236 x3*x17 + 216 x3*x18
      + 108 x3*x19 - 208 x3*x20 - 104 x3*x21 - 32 x3*x22 - 16 x3*x23 + 424 x4^2
      + 424 x4*x5 - 64 x4*x6 - 32 x4*x7 + 400 x4*x8 + 200 x4*x9 + 432 x4*x10
      + 216 x4*x11 + 512 x4*x12 + 256 x4*x13 + 448 x4*x14 + 224 x4*x15
      + 672 x4*x16 + 336 x4*x17 + 432 x4*x18 + 216 x4*x19 + 320 x4*x20
      + 160 x4*x21 + 304 x4*x22 + 152 x4*x23 + 106 x5^2 - 32 x5*x6 - 16 x5*x7
      + 200 x5*x8 + 100 x5*x9 + 216 x5*x10 + 108 x5*x11 + 256 x5*x12
      + 128 x5*x13 + 224 x5*x14 + 112 x5*x15 + 336 x5*x16 + 168 x5*x17
      + 216 x5*x18 + 108 x5*x19 + 160 x5*x20 + 80 x5*x21 + 152 x5*x22
      + 76 x5*x23 + 1304 x6^2 + 1304 x6*x7 + 448 x6*x8 + 224 x6*x9 + 704 x6*x10
      + 352 x6*x11 + 672 x6*x12 + 336 x6*x13 + 560 x6*x14 + 280 x6*x15
      + 480 x6*x16 + 240 x6*x17 + 704 x6*x18 + 352 x6*x19 + 80 x6*x20
      + 40 x6*x21 - 96 x6*x22 - 48 x6*x23 + 326 x7^2 + 224 x7*x8 + 112 x7*x9
      + 352 x7*x10 + 176 x7*x11 + 336 x7*x12 + 168 x7*x13 + 280 x7*x14
      + 140 x7*x15 + 240 x7*x16 + 120 x7*x17 + 352 x7*x18 + 176 x7*x19
      + 40 x7*x20 + 20 x7*x21 - 48 x7*x22 - 24 x7*x23 + 440 x8^2 + 440 x8*x9
      + 320 x8*x10 + 160 x8*x11 + 512 x8*x12 + 256 x8*x13 + 528 x8*x14
      + 264 x8*x15 + 416 x8*x16 + 208 x8*x17 + 400 x8*x18 + 200 x8*x19
      + 336 x8*x20 + 168 x8*x21 + 64 x8*x22 + 32 x8*x23 + 110 x9^2 + 160 x9*x10
      + 80 x9*x11 + 256 x9*x12 + 128 x9*x13 + 264 x9*x14 + 132 x9*x15
      + 208 x9*x16 + 104 x9*x17 + 200 x9*x18 + 100 x9*x19 + 168 x9*x20
      + 84 x9*x21 + 32 x9*x22 + 16 x9*x23 + 504 x10^2 + 504 x10*x11
      + 752 x10*x12 + 376 x10*x13 + 528 x10*x14 + 264 x10*x15 + 864 x10*x16
      + 432 x10*x17 + 272 x10*x18 + 136 x10*x19 - 128 x10*x20 - 64 x10*x21
      + 112 x10*x22 + 56 x10*x23 + 126 x11^2 + 376 x11*x12 + 188 x11*x13
      + 264 x11*x14 + 132 x11*x15 + 432 x11*x16 + 216 x11*x17 + 136 x11*x18
      + 68 x11*x19 - 64 x11*x20 - 32 x11*x21 + 56 x11*x22 + 28 x11*x23
      + 536 x12^2 + 536 x12*x13 + 736 x12*x14 + 368 x12*x15 + 896 x12*x16
      + 448 x12*x17 + 272 x12*x18 + 136 x12*x19 + 224 x12*x20 + 112 x12*x21
      + 64 x12*x22 + 32 x12*x23 + 134 x13^2 + 368 x13*x14 + 184 x13*x15
      + 448 x13*x16 + 224 x13*x17 + 136 x13*x18 + 68 x13*x19 + 112 x13*x20
      + 56 x13*x21 + 32 x13*x22 + 16 x13*x23 + 552 x14^2 + 552 x14*x15
      + 640 x14*x16 + 320 x14*x17 + 272 x14*x18 + 136 x14*x19 + 32 x14*x20
      + 16 x14*x21 - 80 x14*x22 - 40 x14*x23 + 138 x15^2 + 320 x15*x16
      + 160 x15*x17 + 136 x15*x18 + 68 x15*x19 + 16 x15*x20 + 8 x15*x21
      - 40 x15*x22 - 20 x15*x23 + 712 x16^2 + 712 x16*x17 + 320 x16*x18
      + 160 x16*x19 + 336 x16*x20 + 168 x16*x21 + 768 x16*x22 + 384 x16*x23
      + 178 x17^2 + 160 x17*x18 + 80 x17*x19 + 168 x17*x20 + 84 x17*x21
      + 384 x17*x22 + 192 x17*x23 + 608 x18^2 + 608 x18*x19 + 32 x18*x20
      + 16 x18*x21 + 192 x18*x22 + 96 x18*x23 + 152 x19^2 + 16 x19*x20
      + 8 x19*x21 + 96 x19*x22 + 48 x19*x23 + 808 x20^2 + 808 x20*x21
      + 624 x20*x22 + 312 x20*x23 + 202 x21^2 + 312 x21*x22 + 156 x21*x23
      + 1000 x22^2 + 1000 x22*x23 + 250 x23^2 ]/2
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
