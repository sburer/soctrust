Q = [ 2     3    12
      3   -19     6
     12     6     0 ]

c = [ 14
      14
       9 ]

l = -0.5
u =  0.0
eps = 1.2

n = size(Q, 1);

X = sdpvar(n,n,'symm');
x = sdpvar(n,1);
Y = [1,x';x,X];

obj = Q(:)'*X(:) + c'*x;

con = [ Y >= 0 ; trace(X) <= 1 ];
con = [ con ; X(1,1) + l*u + eps*X(2,1) <= (l+u)*x(1) + eps*u*x(2) ];
con = [ con ; norm( X(:,1) - l*x + eps*X(:,2) ) <= x(1) - l + eps*x(2) ];
con = [ con ; norm( u*x - X(:,1) ) <= u - x(1) ];

nost = solvesdp(con,obj);

if nost.problem == 0
    fprintf('optimal value = %.12f\n', double(obj));
    eigvals = sort(eig(double(Y)))
    eigvalratio  = eigvals(n+1)/abs(eigvals(n))
end
