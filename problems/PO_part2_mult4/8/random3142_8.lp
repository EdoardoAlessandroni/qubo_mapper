\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: 230 x0 + 115 x1 + 686 x2 + 343 x3 - 96 x4 - 48 x5 + 628 x6 + 314 x7 + [
      1000 x0^2 + 1000 x0*x1 + 1024 x0*x2 + 512 x0*x3 + 336 x0*x4 + 168 x0*x5
      + 544 x0*x6 + 272 x0*x7 + 250 x1^2 + 512 x1*x2 + 256 x1*x3 + 168 x1*x4
      + 84 x1*x5 + 272 x1*x6 + 136 x1*x7 + 1352 x2^2 + 1352 x2*x3 + 736 x2*x4
      + 368 x2*x5 + 1024 x2*x6 + 512 x2*x7 + 338 x3^2 + 368 x3*x4 + 184 x3*x5
      + 512 x3*x6 + 256 x3*x7 + 512 x4^2 + 512 x4*x5 + 1104 x4*x6 + 552 x4*x7
      + 128 x5^2 + 552 x5*x6 + 276 x5*x7 + 2112 x6^2 + 2112 x6*x7 + 528 x7^2 ]/2
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
