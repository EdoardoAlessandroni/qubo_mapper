\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: 208 x0 + 104 x1 - 856 x2 - 428 x3 + [ 680 x0^2 + 680 x0*x1 + 1824 x0*x2
      + 912 x0*x3 + 170 x1^2 + 912 x1*x2 + 456 x1*x3 + 3232 x2^2 + 3232 x2*x3
      + 808 x3^2 ]/2
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
