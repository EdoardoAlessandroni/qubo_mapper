\ This file has been generated by DOcplex
\ ENCODING=ISO-8859-1
\Problem name: random_problem.lp

Minimize
 obj: - 732 x0 - 366 x1 - 474 x2 - 237 x3 + 510 x4 + 255 x5 + 230 x6 + 115 x7
      + 70 x8 + 35 x9 + [ 632 x0^2 + 632 x0*x1 + 1008 x0*x2 + 504 x0*x3
      + 928 x0*x4 + 464 x0*x5 + 688 x0*x6 + 344 x0*x7 + 224 x0*x8 + 112 x0*x9
      + 158 x1^2 + 504 x1*x2 + 252 x1*x3 + 464 x1*x4 + 232 x1*x5 + 344 x1*x6
      + 172 x1*x7 + 112 x1*x8 + 56 x1*x9 + 1096 x2^2 + 1096 x2*x3 + 1040 x2*x4
      + 520 x2*x5 + 1056 x2*x6 + 528 x2*x7 + 448 x2*x8 + 224 x2*x9 + 274 x3^2
      + 520 x3*x4 + 260 x3*x5 + 528 x3*x6 + 264 x3*x7 + 224 x3*x8 + 112 x3*x9
      + 1064 x4^2 + 1064 x4*x5 + 672 x4*x6 + 336 x4*x7 + 448 x4*x8 + 224 x4*x9
      + 266 x5^2 + 336 x5*x6 + 168 x5*x7 + 224 x5*x8 + 112 x5*x9 + 1000 x6^2
      + 1000 x6*x7 + 640 x6*x8 + 320 x6*x9 + 250 x7^2 + 320 x7*x8 + 160 x7*x9
      + 440 x8^2 + 440 x8*x9 + 110 x9^2 ]/2
Subject To
 c0: 2 x0 + x1 + 2 x2 + x3 + 2 x4 + x5 + 2 x6 + x7 + 2 x8 + x9 = 3

Bounds
 0 <= x0 <= 1
 0 <= x1 <= 1
 0 <= x2 <= 1
 0 <= x3 <= 1
 0 <= x4 <= 1
 0 <= x5 <= 1
 0 <= x6 <= 1
 0 <= x7 <= 1
 0 <= x8 <= 1
 0 <= x9 <= 1

Binaries
 x0 x1 x2 x3 x4 x5 x6 x7 x8 x9
End
