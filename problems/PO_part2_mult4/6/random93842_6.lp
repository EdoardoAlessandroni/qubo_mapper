\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: 686 x0 + 343 x1 - 422 x2 - 211 x3 - 672 x4 - 336 x5 + [ 1352 x0^2
      + 1352 x0*x1 + 1056 x0*x2 + 528 x0*x3 + 800 x0*x4 + 400 x0*x5 + 338 x1^2
      + 528 x1*x2 + 264 x1*x3 + 400 x1*x4 + 200 x1*x5 + 576 x2^2 + 576 x2*x3
      + 448 x2*x4 + 224 x2*x5 + 144 x3^2 + 224 x3*x4 + 112 x3*x5 + 1048 x4^2
      + 1048 x4*x5 + 262 x5^2 ]/2
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
