\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: - 434 x0 - 217 x1 + 374 x2 + 187 x3 - 224 x4 - 112 x5 - 108 x6 - 54 x7 + [
      3328 x0^2 + 3328 x0*x1 + 512 x0*x2 + 256 x0*x3 + 592 x0*x6 + 296 x0*x7
      + 832 x1^2 + 256 x1*x2 + 128 x1*x3 + 296 x1*x6 + 148 x1*x7 + 600 x2^2
      + 600 x2*x3 + 928 x2*x4 + 464 x2*x5 + 656 x2*x6 + 328 x2*x7 + 150 x3^2
      + 464 x3*x4 + 232 x3*x5 + 328 x3*x6 + 164 x3*x7 + 1152 x4^2 + 1152 x4*x5
      + 576 x4*x6 + 288 x4*x7 + 288 x5^2 + 288 x5*x6 + 144 x5*x7 + 360 x6^2
      + 360 x6*x7 + 90 x7^2 ]/2
Subject To
 c0: 2 x0 + x1 + 2 x2 + x3 + 2 x4 + x5 + 2 x6 + x7 = 3

Bounds
 0 <= x0 <= 1
 0 <= x1 <= 1
 0 <= x2 <= 1
 0 <= x3 <= 1
 0 <= x4 <= 1
 0 <= x5 <= 1
 0 <= x6 <= 1
 0 <= x7 <= 1

Binaries
 x0 x1 x2 x3 x4 x5 x6 x7
End
