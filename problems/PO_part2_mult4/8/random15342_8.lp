\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: - 88 x0 - 44 x1 + 200 x2 + 100 x3 - 86 x4 - 43 x5 - 298 x6 - 149 x7 + [
      944 x0^2 + 944 x0*x1 + 896 x0*x2 + 448 x0*x3 + 784 x0*x4 + 392 x0*x5
      + 656 x0*x6 + 328 x0*x7 + 236 x1^2 + 448 x1*x2 + 224 x1*x3 + 392 x1*x4
      + 196 x1*x5 + 328 x1*x6 + 164 x1*x7 + 808 x2^2 + 808 x2*x3 + 720 x2*x4
      + 360 x2*x5 + 880 x2*x6 + 440 x2*x7 + 202 x3^2 + 360 x3*x4 + 180 x3*x5
      + 440 x3*x6 + 220 x3*x7 + 496 x4^2 + 496 x4*x5 + 656 x4*x6 + 328 x4*x7
      + 124 x5^2 + 328 x5*x6 + 164 x5*x7 + 768 x6^2 + 768 x6*x7 + 192 x7^2 ]/2
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
