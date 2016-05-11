reset;
model Section_5_2_example.mod;
data  Section_5_2_example.dat;
option solver couenne;
solve;

display obj;

printf "x(1) = %.15f\n", x[1];
printf "x(2) = %.15f\n", x[2];

printf "objective value = %.12f\n", obj;
