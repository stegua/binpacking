param n;
set I := 1..n;

param m;
set B := 1..m;

param k;
set K := 1..k;

param A {I,K};
param b {K};

param c {I,B} default 0;

var x {I,B} binary;
var s {B} >= 0;

minimize Cost:
   sum {i in I, j in B} c[i,j]*x[i,j];

s.t. Assign {i in I}:
   sum {j in B} x[i,j] = 1;

s.t. Dimensions {j in B, l in K}:
   sum {i in I} A[i,l]*x[i,j] <= b[l] + s[j];

s.t. NoSlack:
   sum {j in B} s[j] = 0;
