\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: - 288 x0 - 144 x1 - 244 x2 - 122 x3 - 106 x4 - 53 x5 - 236 x6 - 118 x7
      + 66 x8 + 33 x9 + [ 320 x0^2 + 320 x0*x1 - 144 x0*x2 - 72 x0*x3
      + 208 x0*x4 + 104 x0*x5 + 112 x0*x6 + 56 x0*x7 - 96 x0*x8 - 48 x0*x9
      + 80 x1^2 - 72 x1*x2 - 36 x1*x3 + 104 x1*x4 + 52 x1*x5 + 56 x1*x6
      + 28 x1*x7 - 48 x1*x8 - 24 x1*x9 + 736 x2^2 + 736 x2*x3 + 448 x2*x4
      + 224 x2*x5 + 64 x2*x6 + 32 x2*x7 + 1024 x2*x8 + 512 x2*x9 + 184 x3^2
      + 224 x3*x4 + 112 x3*x5 + 32 x3*x6 + 16 x3*x7 + 512 x3*x8 + 256 x3*x9
      + 344 x4^2 + 344 x4*x5 + 256 x4*x6 + 128 x4*x7 + 176 x4*x8 + 88 x4*x9
      + 86 x5^2 + 128 x5*x6 + 64 x5*x7 + 88 x5*x8 + 44 x5*x9 + 464 x6^2
      + 464 x6*x7 + 176 x6*x8 + 88 x6*x9 + 116 x7^2 + 88 x7*x8 + 44 x7*x9
      + 944 x8^2 + 944 x8*x9 + 236 x9^2 ]/2
Subject To
 c0: 2 x0 + x1 + 2 x2 + x3 + 2 x4 + x5 + 2 x6 + x7 + 2 x8 + x9 = 3

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

Binaries
 x0 x1 x2 x3 x4 x5 x6 x7 x8 x9
End
