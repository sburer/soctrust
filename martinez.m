function [success,Q,c,r1,H,h,r2,cuts,eigY,LB,UB,sec] = martinez(n,Q,c,r1,H,h,r2)

tic

tol = 1.0e-6;
maxiter = 25;

if nargin == 0

  error('Need at least one argument');

elseif nargin < 7

  %% Set basic radius to be the dimension (Martinez suggests that a
  %% larger radius increases the chances of a local-nonglobal minimum)

  r1 = n;

  %% Set up diagonal Q which is (likely) not psd and for which the
  %% bottom eigenvalue is (likely) unique

  Q = diag(sort(2*rand(n,1)-1));

  %% Setup random c (nothing special here)

  c = 2*rand(n,1) - 1;

  %% Solve regular TRS (only considered successful if SDP delivers
  %% rank-1 solution, which is on the boundary)

  x = sdpvar(n,1);
  X = sdpvar(n,n,'symm');
  Y = [1,x';x,X];
  cons = [ Y >= 0 ; trace(X) <= r1^2 ];
  obj = Q(:)'*X(:) + c'*x;
  solvesdp(cons,obj);

  if rank(double(Y),tol) > 1 | abs(norm(double(x)) - r1) > tol
    success = 0;
    H = [];
    h = [];
    r2 = [];
    cuts = -1;
    eigY = [];
    LB = [];
    UB = [];
    sec = [];
    return
  end

  %% Rotate problem data so that the solution to TRS is actually r1*e1

  [tmp,~] = qr(double(x)); tmp = tmp';
  Q = tmp'*Q*tmp;
  c = tmp'*c;

  %% Setup a random coordinate-aligned ellipsoid, which does not include
  %% r1*e1 (i.e., it cuts off the global solution of TRS)

  r2 = r1;
  H = eye(n);
  H(1,1) = 2; % This ensures r1*e1 is not a solution
  for i = 2:n
    H(i,i) = 0.5 + (2-0.5)*rand; % Random between 0.5 and 2
  end
%   [V,D] = eig(H); d = diag(D); Hsqrt = V*diag(sqrt(d))*V';
  h = zeros(n,1);
%   h = (2*rand(n,1)-1)/n; % Relatively small shift of the center. Could this be a problem?

end

%% Check that n matches the size of Q (this is not a very thorough input
%% check)

if n ~= size(Q,1)
  error('Size mismatch');
end

%% Derive the sqrt of H

[V,D] = eig(H); d = diag(D); Hsqrt = V*diag(sqrt(d))*V';

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

cuts = 0;

for iter = 1:maxiter

  %% Solve current relaxation

  solvesdp(cons,obj);
  fprintf('iter = %d   relaxed value = %f\n', iter, double(obj));
%   pause

  %% Store portion of solution

  dx = double(x);
  dX = double(X);
  eigY = eig(double(Y));

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
%   double(sepobj);

  %% If no separation, quit

  if double(sepobj)/(norm(Q,'fro') + norm(c)) > -tol
%     fprintf('Solved!\n');
%     eig(double(Y))
%     double(x)
    break
  end


  %% Add SOC-RLT constraint to relaxation

  da = double(a);
%   pause

  cons = [ cons ; norm( Hsqrt*(r1*(x-h) - X*da + h*da'*x )) <= r2*(r1 - da'*x) ];

  cuts = cuts+1;

end

LB = double(obj);
UB = dx'*Q*dx + c'*dx; 
success = 1;
sec = toc;
