\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: 40 x0 + 20 x1 - 300 x2 - 150 x3 - 468 x4 - 234 x5 - 152 x6 - 76 x7
      + 370 x8 + 185 x9 + 364 x10 + 182 x11 - 698 x12 - 349 x13 - 78 x14
      - 39 x15 + 376 x16 + 188 x17 - 242 x18 - 121 x19 + 10 x20 + 5 x21
      + 104 x22 + 52 x23 + [ 984 x0^2 + 984 x0*x1 + 1760 x0*x2 + 880 x0*x3
      + 464 x0*x4 + 232 x0*x5 + 1040 x0*x6 + 520 x0*x7 + 880 x0*x8 + 440 x0*x9
      + 432 x0*x10 + 216 x0*x11 + 352 x0*x12 + 176 x0*x13 + 864 x0*x14
      + 432 x0*x15 + 992 x0*x16 + 496 x0*x17 + 848 x0*x18 + 424 x0*x19
      + 480 x0*x20 + 240 x0*x21 + 1248 x0*x22 + 624 x0*x23 + 246 x1^2
      + 880 x1*x2 + 440 x1*x3 + 232 x1*x4 + 116 x1*x5 + 520 x1*x6 + 260 x1*x7
      + 440 x1*x8 + 220 x1*x9 + 216 x1*x10 + 108 x1*x11 + 176 x1*x12 + 88 x1*x13
      + 432 x1*x14 + 216 x1*x15 + 496 x1*x16 + 248 x1*x17 + 424 x1*x18
      + 212 x1*x19 + 240 x1*x20 + 120 x1*x21 + 624 x1*x22 + 312 x1*x23
      + 1600 x2^2 + 1600 x2*x3 + 32 x2*x4 + 16 x2*x5 + 784 x2*x6 + 392 x2*x7
      + 544 x2*x8 + 272 x2*x9 + 656 x2*x10 + 328 x2*x11 + 688 x2*x14
      + 344 x2*x15 + 784 x2*x16 + 392 x2*x17 + 880 x2*x18 + 440 x2*x19
      + 128 x2*x20 + 64 x2*x21 + 320 x2*x22 + 160 x2*x23 + 400 x3^2 + 16 x3*x4
      + 8 x3*x5 + 392 x3*x6 + 196 x3*x7 + 272 x3*x8 + 136 x3*x9 + 328 x3*x10
      + 164 x3*x11 + 344 x3*x14 + 172 x3*x15 + 392 x3*x16 + 196 x3*x17
      + 440 x3*x18 + 220 x3*x19 + 64 x3*x20 + 32 x3*x21 + 160 x3*x22 + 80 x3*x23
      + 504 x4^2 + 504 x4*x5 + 272 x4*x6 + 136 x4*x7 + 784 x4*x8 + 392 x4*x9
      - 32 x4*x10 - 16 x4*x11 + 496 x4*x12 + 248 x4*x13 + 432 x4*x14
      + 216 x4*x15 - 48 x4*x16 - 24 x4*x17 + 464 x4*x18 + 232 x4*x19
      + 768 x4*x20 + 384 x4*x21 + 1552 x4*x22 + 776 x4*x23 + 126 x5^2
      + 136 x5*x6 + 68 x5*x7 + 392 x5*x8 + 196 x5*x9 - 16 x5*x10 - 8 x5*x11
      + 248 x5*x12 + 124 x5*x13 + 216 x5*x14 + 108 x5*x15 - 24 x5*x16
      - 12 x5*x17 + 232 x5*x18 + 116 x5*x19 + 384 x5*x20 + 192 x5*x21
      + 776 x5*x22 + 388 x5*x23 + 480 x6^2 + 480 x6*x7 + 720 x6*x8 + 360 x6*x9
      + 272 x6*x10 + 136 x6*x11 + 160 x6*x12 + 80 x6*x13 + 384 x6*x14
      + 192 x6*x15 + 816 x6*x16 + 408 x6*x17 + 464 x6*x18 + 232 x6*x19
      + 128 x6*x20 + 64 x6*x21 + 544 x6*x22 + 272 x6*x23 + 120 x7^2 + 360 x7*x8
      + 180 x7*x9 + 136 x7*x10 + 68 x7*x11 + 80 x7*x12 + 40 x7*x13 + 192 x7*x14
      + 96 x7*x15 + 408 x7*x16 + 204 x7*x17 + 232 x7*x18 + 116 x7*x19
      + 64 x7*x20 + 32 x7*x21 + 272 x7*x22 + 136 x7*x23 + 960 x8^2 + 960 x8*x9
      + 144 x8*x10 + 72 x8*x11 + 416 x8*x12 + 208 x8*x13 + 624 x8*x14
      + 312 x8*x15 + 416 x8*x16 + 208 x8*x17 + 896 x8*x18 + 448 x8*x19
      + 480 x8*x20 + 240 x8*x21 + 1072 x8*x22 + 536 x8*x23 + 240 x9^2
      + 72 x9*x10 + 36 x9*x11 + 208 x9*x12 + 104 x9*x13 + 312 x9*x14
      + 156 x9*x15 + 208 x9*x16 + 104 x9*x17 + 448 x9*x18 + 224 x9*x19
      + 240 x9*x20 + 120 x9*x21 + 536 x9*x22 + 268 x9*x23 + 608 x10^2
      + 608 x10*x11 + 336 x10*x14 + 168 x10*x15 + 496 x10*x16 + 248 x10*x17
      + 304 x10*x18 + 152 x10*x19 + 192 x10*x20 + 96 x10*x21 + 256 x10*x22
      + 128 x10*x23 + 152 x11^2 + 168 x11*x14 + 84 x11*x15 + 248 x11*x16
      + 124 x11*x17 + 152 x11*x18 + 76 x11*x19 + 96 x11*x20 + 48 x11*x21
      + 128 x11*x22 + 64 x11*x23 + 472 x12^2 + 472 x12*x13 + 256 x12*x14
      + 128 x12*x15 + 64 x12*x16 + 32 x12*x17 + 528 x12*x18 + 264 x12*x19
      + 192 x12*x20 + 96 x12*x21 + 560 x12*x22 + 280 x12*x23 + 118 x13^2
      + 128 x13*x14 + 64 x13*x15 + 32 x13*x16 + 16 x13*x17 + 264 x13*x18
      + 132 x13*x19 + 96 x13*x20 + 48 x13*x21 + 280 x13*x22 + 140 x13*x23
      + 432 x14^2 + 432 x14*x15 + 512 x14*x16 + 256 x14*x17 + 704 x14*x18
      + 352 x14*x19 + 768 x14*x20 + 384 x14*x21 + 1120 x14*x22 + 560 x14*x23
      + 108 x15^2 + 256 x15*x16 + 128 x15*x17 + 352 x15*x18 + 176 x15*x19
      + 384 x15*x20 + 192 x15*x21 + 560 x15*x22 + 280 x15*x23 + 792 x16^2
      + 792 x16*x17 + 320 x16*x18 + 160 x16*x19 + 16 x16*x20 + 8 x16*x21
      + 368 x16*x22 + 184 x16*x23 + 198 x17^2 + 160 x17*x18 + 80 x17*x19
      + 8 x17*x20 + 4 x17*x21 + 184 x17*x22 + 92 x17*x23 + 672 x18^2
      + 672 x18*x19 + 608 x18*x20 + 304 x18*x21 + 752 x18*x22 + 376 x18*x23
      + 168 x19^2 + 304 x19*x20 + 152 x19*x21 + 376 x19*x22 + 188 x19*x23
      + 1000 x20^2 + 1000 x20*x21 + 2512 x20*x22 + 1256 x20*x23 + 250 x21^2
      + 1256 x21*x22 + 628 x21*x23 + 2416 x22^2 + 2416 x22*x23 + 604 x23^2 ]/2
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
