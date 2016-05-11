reset;
model TTRS.mod;
data  TTRS.dat;
option solver couenne;
solve;

# display obj;
# printf "x(1) = %.15f\n", x[1];
# printf "x(2) = %.15f\n", x[2];
# printf "x(3) = %.15f\n", x[3];

# printf "keep %.12f\n", obj;
