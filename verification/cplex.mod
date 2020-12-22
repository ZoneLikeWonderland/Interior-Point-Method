param m ;
param n ;
set N := 1..n;
set M := 1..m;
param A {M, N};
param c{N};
param b {M};


var x {N} >=0;
 
minimize Cost: sum {i in N} c[i]*x[i];

subject to C{i in M}: sum {j in N } A[i,j]*x[j]=b[i];