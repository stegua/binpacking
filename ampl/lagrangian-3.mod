param n;
set I := 1..n;

param m;
set B := 1..m;

param k;
set K := 1..k;

param A {I,K};
param b {K};

param c {I,B} default 0;

param alpha {i in I, j in B} default 1;
param G  {i in I, j in B} default 0;
param H1 {i in I, j in B} default 0;
param H2 {i in I, j in B} default 0;
param f default 20.0;
param t default 0.0;

var x {I,B} binary;
var y {I,B} binary;

# Rilassamento lineare
minimize P1:
   sum {i in I, j in B} c[i,j]*x[i,j] - sum {i in I, j in B} alpha[i,j]*x[i,j];

s.t. Assign {i in I}:
   sum {j in B} x[i,j] = 1;

minimize P2:
   sum {i in I, j in B} alpha[i,j]*y[i,j];

s.t. Dimensions {j in B, l in K}:
   sum {i in I} A[i,l]*y[i,j] <= b[l];


