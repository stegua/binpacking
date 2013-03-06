param n;
set I := 1..n;

param m;
set B := 1..m;

param k;
set K := 1..k;

param A {I,K};
param b {K};

param c {I,B} default 0;

param alpha {j in B, l in K} >= 0 default 0;
param G  {j in B, l in K} default 0;
param H1 {j in B, l in K} default 0;
param H2 {j in B, l in K} default 0;
param f default 20.0;
param t default 0.0;

var x {I,B} >= 0, <= 1;

minimize cLB:
   sum {i in I, j in B} c[i,j]*x[i,j] - sum {j in B, l in K} alpha[j,l]*(b[l] - sum {i in I} A[i,l]*x[i,j]);

s.t. Assign {i in I}:
   sum {j in B} x[i,j] = 1;

