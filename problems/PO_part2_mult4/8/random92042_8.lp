\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: 52 x0 + 26 x1 - 422 x2 - 211 x3 + 18 x4 + 9 x5 - 310 x6 - 155 x7 + [
      2448 x0^2 + 2448 x0*x1 + 1792 x0*x2 + 896 x0*x3 + 1888 x0*x4 + 944 x0*x5
      + 928 x0*x6 + 464 x0*x7 + 612 x1^2 + 896 x1*x2 + 448 x1*x3 + 944 x1*x4
      + 472 x1*x5 + 464 x1*x6 + 232 x1*x7 + 576 x2^2 + 576 x2*x3 + 912 x2*x4
      + 456 x2*x5 + 784 x2*x6 + 392 x2*x7 + 144 x3^2 + 456 x3*x4 + 228 x3*x5
      + 392 x3*x6 + 196 x3*x7 + 832 x4^2 + 832 x4*x5 + 496 x4*x6 + 248 x4*x7
      + 208 x5^2 + 248 x5*x6 + 124 x5*x7 + 664 x6^2 + 664 x6*x7 + 166 x7^2 ]/2
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
