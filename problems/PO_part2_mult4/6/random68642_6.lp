\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: - 68 x0 - 34 x1 + 36 x2 + 18 x3 - 2 x4 - x5 + [ 2160 x0^2 + 2160 x0*x1
      + 1408 x0*x2 + 704 x0*x3 + 1408 x0*x4 + 704 x0*x5 + 540 x1^2 + 704 x1*x2
      + 352 x1*x3 + 704 x1*x4 + 352 x1*x5 + 712 x2^2 + 712 x2*x3 + 1040 x2*x4
      + 520 x2*x5 + 178 x3^2 + 520 x3*x4 + 260 x3*x5 + 2128 x4^2 + 2128 x4*x5
      + 532 x5^2 ]/2
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
