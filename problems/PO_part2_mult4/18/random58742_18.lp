\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: - 2 x0 - x1 - 104 x2 - 52 x3 - 244 x4 - 122 x5 - 474 x6 - 237 x7 + 328 x8
      + 164 x9 + 628 x10 + 314 x11 - 154 x12 - 77 x13 - 88 x14 - 44 x15
      + 498 x16 + 249 x17 + [ 2128 x0^2 + 2128 x0*x1 - 512 x0*x2 - 256 x0*x3
      + 368 x0*x4 + 184 x0*x5 + 1088 x0*x6 + 544 x0*x7 + 1152 x0*x8 + 576 x0*x9
      + 528 x0*x10 + 264 x0*x11 - 720 x0*x12 - 360 x0*x13 + 416 x0*x14
      + 208 x0*x15 + 432 x0*x16 + 216 x0*x17 + 532 x1^2 - 256 x1*x2 - 128 x1*x3
      + 184 x1*x4 + 92 x1*x5 + 544 x1*x6 + 272 x1*x7 + 576 x1*x8 + 288 x1*x9
      + 264 x1*x10 + 132 x1*x11 - 360 x1*x12 - 180 x1*x13 + 208 x1*x14
      + 104 x1*x15 + 216 x1*x16 + 108 x1*x17 + 2608 x2^2 + 2608 x2*x3
      + 1600 x2*x4 + 800 x2*x5 + 2288 x2*x6 + 1144 x2*x7 + 304 x2*x8 + 152 x2*x9
      + 1888 x2*x10 + 944 x2*x11 + 3184 x2*x12 + 1592 x2*x13 + 1152 x2*x14
      + 576 x2*x15 + 1712 x2*x16 + 856 x2*x17 + 652 x3^2 + 800 x3*x4 + 400 x3*x5
      + 1144 x3*x6 + 572 x3*x7 + 152 x3*x8 + 76 x3*x9 + 944 x3*x10 + 472 x3*x11
      + 1592 x3*x12 + 796 x3*x13 + 576 x3*x14 + 288 x3*x15 + 856 x3*x16
      + 428 x3*x17 + 736 x4^2 + 736 x4*x5 + 1312 x4*x6 + 656 x4*x7 + 832 x4*x8
      + 416 x4*x9 + 912 x4*x10 + 456 x4*x11 + 976 x4*x12 + 488 x4*x13
      + 672 x4*x14 + 336 x4*x15 + 992 x4*x16 + 496 x4*x17 + 184 x5^2 + 656 x5*x6
      + 328 x5*x7 + 416 x5*x8 + 208 x5*x9 + 456 x5*x10 + 228 x5*x11 + 488 x5*x12
      + 244 x5*x13 + 336 x5*x14 + 168 x5*x15 + 496 x5*x16 + 248 x5*x17
      + 1096 x6^2 + 1096 x6*x7 + 496 x6*x8 + 248 x6*x9 + 1424 x6*x10
      + 712 x6*x11 + 1520 x6*x12 + 760 x6*x13 + 880 x6*x14 + 440 x6*x15
      + 1088 x6*x16 + 544 x6*x17 + 274 x7^2 + 248 x7*x8 + 124 x7*x9 + 712 x7*x10
      + 356 x7*x11 + 760 x7*x12 + 380 x7*x13 + 440 x7*x14 + 220 x7*x15
      + 544 x7*x16 + 272 x7*x17 + 1128 x8^2 + 1128 x8*x9 + 688 x8*x10
      + 344 x8*x11 - 64 x8*x12 - 32 x8*x13 - 128 x8*x14 - 64 x8*x15 + 432 x8*x16
      + 216 x8*x17 + 282 x9^2 + 344 x9*x10 + 172 x9*x11 - 32 x9*x12 - 16 x9*x13
      - 64 x9*x14 - 32 x9*x15 + 216 x9*x16 + 108 x9*x17 + 2112 x10^2
      + 2112 x10*x11 + 1104 x10*x12 + 552 x10*x13 + 208 x10*x14 + 104 x10*x15
      + 1184 x10*x16 + 592 x10*x17 + 528 x11^2 + 552 x11*x12 + 276 x11*x13
      + 104 x11*x14 + 52 x11*x15 + 592 x11*x16 + 296 x11*x17 + 1264 x12^2
      + 1264 x12*x13 + 768 x12*x14 + 384 x12*x15 + 960 x12*x16 + 480 x12*x17
      + 316 x13^2 + 384 x13*x14 + 192 x13*x15 + 480 x13*x16 + 240 x13*x17
      + 472 x14^2 + 472 x14*x15 + 672 x14*x16 + 336 x14*x17 + 118 x15^2
      + 336 x15*x16 + 168 x15*x17 + 680 x16^2 + 680 x16*x17 + 170 x17^2 ]/2
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
