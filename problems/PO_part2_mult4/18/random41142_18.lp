\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: - 198 x0 - 99 x1 - 156 x2 - 78 x3 + 70 x4 + 35 x5 + 370 x6 + 185 x7
      + 200 x8 + 100 x9 - 152 x10 - 76 x11 - 258 x12 - 129 x13 - 1014 x14
      - 507 x15 + 258 x16 + 129 x17 + [ 464 x0^2 + 464 x0*x1 + 336 x0*x2
      + 168 x0*x3 + 400 x0*x4 + 200 x0*x5 + 640 x0*x6 + 320 x0*x7 + 336 x0*x8
      + 168 x0*x9 + 368 x0*x10 + 184 x0*x11 - 48 x0*x12 - 24 x0*x13 + 576 x0*x14
      + 288 x0*x15 + 320 x0*x16 + 160 x0*x17 + 116 x1^2 + 168 x1*x2 + 84 x1*x3
      + 200 x1*x4 + 100 x1*x5 + 320 x1*x6 + 160 x1*x7 + 168 x1*x8 + 84 x1*x9
      + 184 x1*x10 + 92 x1*x11 - 24 x1*x12 - 12 x1*x13 + 288 x1*x14 + 144 x1*x15
      + 160 x1*x16 + 80 x1*x17 + 504 x2^2 + 504 x2*x3 + 320 x2*x4 + 160 x2*x5
      + 720 x2*x6 + 360 x2*x7 + 512 x2*x8 + 256 x2*x9 + 976 x2*x10 + 488 x2*x11
      - 128 x2*x12 - 64 x2*x13 + 1504 x2*x14 + 752 x2*x15 + 976 x2*x16
      + 488 x2*x17 + 126 x3^2 + 160 x3*x4 + 80 x3*x5 + 360 x3*x6 + 180 x3*x7
      + 256 x3*x8 + 128 x3*x9 + 488 x3*x10 + 244 x3*x11 - 64 x3*x12 - 32 x3*x13
      + 752 x3*x14 + 376 x3*x15 + 488 x3*x16 + 244 x3*x17 + 440 x4^2 + 440 x4*x5
      + 320 x4*x6 + 160 x4*x7 + 160 x4*x8 + 80 x4*x9 + 304 x4*x10 + 152 x4*x11
      + 336 x4*x12 + 168 x4*x13 + 1248 x4*x14 + 624 x4*x15 + 368 x4*x16
      + 184 x4*x17 + 110 x5^2 + 160 x5*x6 + 80 x5*x7 + 80 x5*x8 + 40 x5*x9
      + 152 x5*x10 + 76 x5*x11 + 168 x5*x12 + 84 x5*x13 + 624 x5*x14
      + 312 x5*x15 + 184 x5*x16 + 92 x5*x17 + 960 x6^2 + 960 x6*x7 + 560 x6*x8
      + 280 x6*x9 + 720 x6*x10 + 360 x6*x11 - 48 x6*x12 - 24 x6*x13
      + 1472 x6*x14 + 736 x6*x15 + 1120 x6*x16 + 560 x6*x17 + 240 x7^2
      + 280 x7*x8 + 140 x7*x9 + 360 x7*x10 + 180 x7*x11 - 24 x7*x12 - 12 x7*x13
      + 736 x7*x14 + 368 x7*x15 + 560 x7*x16 + 280 x7*x17 + 808 x8^2 + 808 x8*x9
      + 496 x8*x10 + 248 x8*x11 + 16 x8*x12 + 8 x8*x13 + 400 x8*x14 + 200 x8*x15
      - 48 x8*x16 - 24 x8*x17 + 202 x9^2 + 248 x9*x10 + 124 x9*x11 + 8 x9*x12
      + 4 x9*x13 + 200 x9*x14 + 100 x9*x15 - 24 x9*x16 - 12 x9*x17 + 480 x10^2
      + 480 x10*x11 - 160 x10*x12 - 80 x10*x13 + 1520 x10*x14 + 760 x10*x15
      + 1008 x10*x16 + 504 x10*x17 + 120 x11^2 - 80 x11*x12 - 40 x11*x13
      + 760 x11*x14 + 380 x11*x15 + 504 x11*x16 + 252 x11*x17 + 808 x12^2
      + 808 x12*x13 - 1168 x12*x14 - 584 x12*x15 + 656 x12*x16 + 328 x12*x17
      + 202 x13^2 - 584 x13*x14 - 292 x13*x15 + 328 x13*x16 + 164 x13*x17
      + 4880 x14^2 + 4880 x14*x15 + 2384 x14*x16 + 1192 x14*x17 + 1220 x15^2
      + 1192 x15*x16 + 596 x15*x17 + 1944 x16^2 + 1944 x16*x17 + 486 x17^2 ]/2
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
