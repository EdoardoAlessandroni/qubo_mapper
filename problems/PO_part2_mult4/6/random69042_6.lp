\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: - 128 x0 - 64 x1 + 36 x2 + 18 x3 + 110 x4 + 55 x5 + [ 368 x0^2 + 368 x0*x1
      + 320 x0*x2 + 160 x0*x3 + 92 x1^2 + 160 x1*x2 + 80 x1*x3 + 712 x2^2
      + 712 x2*x3 + 624 x2*x4 + 312 x2*x5 + 178 x3^2 + 312 x3*x4 + 156 x3*x5
      + 1192 x4^2 + 1192 x4*x5 + 298 x5^2 ]/2
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
