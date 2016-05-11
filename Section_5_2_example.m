n = 2;

Q = [ -4 1 
       1 -2 ];

c = [ 1 1 ]';

H = 0.5*[ 3 0 ; 0 1 ];
[V,D] = eig(H); d = diag(D); Hsqrt = V*diag(sqrt(d))*V';
h = [ 0 0 ]';

r1 = 1;
r2 = 1;

%% Setup first SDP

x = sdpvar(n,1);
X = sdpvar(n,n,'symm');
Y = [1,x';x,X];

cons = [ Y >= 0 ];
cons = [ cons ; trace(X) <= r1^2 ];
cons = [ cons ; H(:)'*X(:) - 2*h'*H*x + h'*H*h <= r2^2 ];

obj = Q(:)'*X(:) + c'*x;

%% Setup separation problem (objective TBD)

a = sdpvar(n,1);
A = sdpvar(n,n,'symm');
sepcons = [ [1,a';a,A] >= 0 ];
sepcons = [ sepcons ; trace(A) <= 1 ];

%% Iterate

for iter = 1:100

    %% Solve current relaxation

    solvesdp(cons,obj);
    fprintf('iter = %d   relaxed value = %f   #SOCRLT = %d\n', iter, double(obj), iter-1);

    double(x)
    double(X)
    eig(double(Y))

    %% Store portion of solution

    dx = double(x);
    dX = double(X);

    %% Setup objective for SOC-RLT separation

    tmp1 = r1*(dx-h);
    tmp2 = h*dx' - dX;

    sepobj1 = (r2^2)*(r1^2 - 2*r1*a'*dx + dx'*A*dx);

    sepobj21 = tmp1'*H*tmp1;
    sepobj22 = tmp1'*H*tmp2;
    sepobj23 = tmp2'*H*tmp2;

    sepobj = sepobj1 - (sepobj21 + 2*sepobj22*a + sepobj23(:)'*A(:));

    %% Solve separation subproblem

    solvesdp(sepcons,sepobj);
    double(sepobj)

    %% If no separation, quit

    if double(sepobj) > -1.0e-8
      fprintf('Solved!\n');
      eig(double(Y))
      double(x)
      break
    end

    %% Add SOC-RLT constraint to relaxation

    da = double(a)

    cons = [ cons ; norm( Hsqrt*(r1*(x-h) - X*da + h*da'*x )) <= r2*(r1 - da'*x) ];

end
    
fprintf('iter = %d   relaxed value = %f   #SOCRLT = %d\n', iter, double(obj), iter-1);
