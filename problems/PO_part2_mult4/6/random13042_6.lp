\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: 628 x0 + 314 x1 - 244 x2 - 122 x3 + 562 x4 + 281 x5 + [ 2112 x0^2
      + 2112 x0*x1 + 912 x0*x2 + 456 x0*x3 + 512 x0*x4 + 256 x0*x5 + 528 x1^2
      + 456 x1*x2 + 228 x1*x3 + 256 x1*x4 + 128 x1*x5 + 736 x2^2 + 736 x2*x3
      + 560 x2*x4 + 280 x2*x5 + 184 x3^2 + 280 x3*x4 + 140 x3*x5 + 608 x4^2
      + 608 x4*x5 + 152 x5^2 ]/2
Subject To
 c0: 2 x0 + x1 + 2 x2 + x3 + 2 x4 + x5 = 3

Bounds
 0 <= x0 <= 1
 0 <= x1 <= 1
 0 <= x2 <= 1
 0 <= x3 <= 1
 0 <= x4 <= 1
 0 <= x5 <= 1

Binaries
 x0 x1 x2 x3 x4 x5
End
