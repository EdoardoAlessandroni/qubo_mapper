\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: - 78 x0 - 39 x1 - 200 x2 - 100 x3 + 66 x4 + 33 x5 - 852 x6 - 426 x7
      - 258 x8 - 129 x9 - 68 x10 - 34 x11 - 244 x12 - 122 x13 + 90 x14 + 45 x15
      - 418 x16 - 209 x17 + [ 976 x0^2 + 976 x0*x1 + 96 x0*x2 + 48 x0*x3
      + 720 x0*x4 + 360 x0*x5 + 1040 x0*x6 + 520 x0*x7 + 656 x0*x8 + 328 x0*x9
      + 1008 x0*x10 + 504 x0*x11 + 720 x0*x12 + 360 x0*x13 + 1056 x0*x14
      + 528 x0*x15 + 992 x0*x16 + 496 x0*x17 + 244 x1^2 + 48 x1*x2 + 24 x1*x3
      + 360 x1*x4 + 180 x1*x5 + 520 x1*x6 + 260 x1*x7 + 328 x1*x8 + 164 x1*x9
      + 504 x1*x10 + 252 x1*x11 + 360 x1*x12 + 180 x1*x13 + 528 x1*x14
      + 264 x1*x15 + 496 x1*x16 + 248 x1*x17 + 248 x2^2 + 248 x2*x3 + 416 x2*x4
      + 208 x2*x5 - 48 x2*x6 - 24 x2*x7 + 240 x2*x8 + 120 x2*x9 - 352 x2*x10
      - 176 x2*x11 + 256 x2*x12 + 128 x2*x13 - 448 x2*x14 - 224 x2*x15
      + 320 x2*x16 + 160 x2*x17 + 62 x3^2 + 208 x3*x4 + 104 x3*x5 - 24 x3*x6
      - 12 x3*x7 + 120 x3*x8 + 60 x3*x9 - 176 x3*x10 - 88 x3*x11 + 128 x3*x12
      + 64 x3*x13 - 224 x3*x14 - 112 x3*x15 + 160 x3*x16 + 80 x3*x17 + 944 x4^2
      + 944 x4*x5 + 688 x4*x6 + 344 x4*x7 + 400 x4*x8 + 200 x4*x9 + 160 x4*x10
      + 80 x4*x11 + 1024 x4*x12 + 512 x4*x13 + 32 x4*x14 + 16 x4*x15
      + 608 x4*x16 + 304 x4*x17 + 236 x5^2 + 344 x5*x6 + 172 x5*x7 + 200 x5*x8
      + 100 x5*x9 + 80 x5*x10 + 40 x5*x11 + 512 x5*x12 + 256 x5*x13 + 16 x5*x14
      + 8 x5*x15 + 304 x5*x16 + 152 x5*x17 + 1320 x6^2 + 1320 x6*x7 + 208 x6*x8
      + 104 x6*x9 + 2560 x6*x10 + 1280 x6*x11 + 1008 x6*x12 + 504 x6*x13
      + 2192 x6*x14 + 1096 x6*x15 - 80 x6*x16 - 40 x6*x17 + 330 x7^2 + 104 x7*x8
      + 52 x7*x9 + 1280 x7*x10 + 640 x7*x11 + 504 x7*x12 + 252 x7*x13
      + 1096 x7*x14 + 548 x7*x15 - 40 x7*x16 - 20 x7*x17 + 808 x8^2 + 808 x8*x9
      + 112 x8*x10 + 56 x8*x11 + 16 x8*x12 + 8 x8*x13 + 480 x8*x14 + 240 x8*x15
      + 640 x8*x16 + 320 x8*x17 + 202 x9^2 + 56 x9*x10 + 28 x9*x11 + 8 x9*x12
      + 4 x9*x13 + 240 x9*x14 + 120 x9*x15 + 320 x9*x16 + 160 x9*x17
      + 2160 x10^2 + 2160 x10*x11 + 1024 x10*x12 + 512 x10*x13 + 3024 x10*x14
      + 1512 x10*x15 + 368 x10*x16 + 184 x10*x17 + 540 x11^2 + 512 x11*x12
      + 256 x11*x13 + 1512 x11*x14 + 756 x11*x15 + 184 x11*x16 + 92 x11*x17
      + 736 x12^2 + 736 x12*x13 + 880 x12*x14 + 440 x12*x15 + 768 x12*x16
      + 384 x12*x17 + 184 x13^2 + 440 x13*x14 + 220 x13*x15 + 384 x13*x16
      + 192 x13*x17 + 2528 x14^2 + 2528 x14*x15 + 1232 x14*x16 + 616 x14*x17
      + 632 x15^2 + 616 x15*x16 + 308 x15*x17 + 2128 x16^2 + 2128 x16*x17
      + 532 x17^2 ]/2
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
