\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: - 290 x0 - 145 x1 - 434 x2 - 217 x3 - 198 x4 - 99 x5 + [ 752 x0^2
      + 752 x0*x1 + 688 x0*x2 + 344 x0*x3 + 400 x0*x4 + 200 x0*x5 + 188 x1^2
      + 344 x1*x2 + 172 x1*x3 + 200 x1*x4 + 100 x1*x5 + 3328 x2^2 + 3328 x2*x3
      - 96 x2*x4 - 48 x2*x5 + 832 x3^2 - 48 x3*x4 - 24 x3*x5 + 464 x4^2
      + 464 x4*x5 + 116 x5^2 ]/2
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
