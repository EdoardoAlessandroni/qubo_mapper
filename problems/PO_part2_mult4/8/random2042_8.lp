\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: 56 x0 + 28 x1 + 90 x2 + 45 x3 + 258 x4 + 129 x5 + 584 x6 + 292 x7 + [
      488 x0^2 + 488 x0*x1 + 480 x0*x2 + 240 x0*x3 + 592 x0*x4 + 296 x0*x5
      + 304 x0*x6 + 152 x0*x7 + 122 x1^2 + 240 x1*x2 + 120 x1*x3 + 296 x1*x4
      + 148 x1*x5 + 152 x1*x6 + 76 x1*x7 + 448 x2^2 + 448 x2*x3 + 768 x2*x4
      + 384 x2*x5 + 528 x2*x6 + 264 x2*x7 + 112 x3^2 + 384 x3*x4 + 192 x3*x5
      + 264 x3*x6 + 132 x3*x7 + 1944 x4^2 + 1944 x4*x5 - 32 x4*x6 - 16 x4*x7
      + 486 x5^2 - 16 x5*x6 - 8 x5*x7 + 1304 x6^2 + 1304 x6*x7 + 326 x7^2 ]/2
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
