\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: 10 x0 + 5 x1 + 400 x2 + 200 x3 - 458 x4 - 229 x5 + [ 1000 x0^2
      + 1000 x0*x1 + 2784 x0*x2 + 1392 x0*x3 + 160 x0*x4 + 80 x0*x5 + 250 x1^2
      + 1392 x1*x2 + 696 x1*x3 + 80 x1*x4 + 40 x1*x5 + 2528 x2^2 + 2528 x2*x3
      - 128 x2*x4 - 64 x2*x5 + 632 x3^2 - 64 x3*x4 - 32 x3*x5 + 816 x4^2
      + 816 x4*x5 + 204 x5^2 ]/2
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
