set J := 1..3;

param Q {J,J};
param c {J};
param l;
param u;
param eps;

var x{J};

minimize obj: sum {j in J} c[j]*x[j] + sum {i in J, j in J} Q[i,j]*x[i]*x[j];

s.t. left:  l <= x[1] + eps*x[2];
s.t. right: x[1] <= u;
s.t. ball:  x[1]^2 + x[2]^2 + x[3]^2 <= 1;
