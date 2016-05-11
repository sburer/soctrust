set J;

param Q {J,J};
param c {J};
param H {J}; # Assumed diagonal!
# param h {J}; Assumed zero!!
param r1;
param r2;

var x{J};

minimize obj: sum {j in J} c[j]*x[j] + sum {i in J, j in J} Q[i,j]*x[i]*x[j];

s.t. ball:       sum {j in J} x[j]^2 <= r1^2;
s.t. ellipsoid:  sum {j in J} H[j]*x[j]^2 <= r2^2;
