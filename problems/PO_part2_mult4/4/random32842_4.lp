\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: 10 x0 + 5 x1 - 852 x2 - 426 x3 + [ 1000 x0^2 + 1000 x0*x1 + 752 x0*x2
      + 376 x0*x3 + 250 x1^2 + 376 x1*x2 + 188 x1*x3 + 1320 x2^2 + 1320 x2*x3
      + 330 x3^2 ]/2
Subject To
 c0: 2 x0 + x1 + 2 x2 + x3 = 3

Bounds
 0 <= x0 <= 1
 0 <= x1 <= 1
 0 <= x2 <= 1
 0 <= x3 <= 1

Binaries
 x0 x1 x2 x3
End
