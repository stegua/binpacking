param n;
set I := 1..n;

param m;
set B := 1..m;

param k;
set K := 1..k;

param A {I,K};
param b {K};

param r symbolic in K;

param c {I,B} default 0;

var x {I,B} >= 0, <= 1;

maximize cUB:
   sum {i in I, j in B} c[i,j]*x[i,j];

minimize cLB:
   sum {i in I} A[i,r]*x[i,1];

s.t. Assign {i in I}:
   sum {j in B} x[i,j] = 1;

s.t. Dimensions {j in B, l in K: not (j = 1 and l = r) }:
   sum {i in I} A[i,l]*x[i,j] <= b[l];

