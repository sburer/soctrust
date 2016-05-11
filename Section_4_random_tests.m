n = 3;

% ratios = [];

minratio = +Inf;

for iter = 1:50000

  u = -1 + rand*2;
  l = -1 + rand*(u+1);

  if l >= u
    error('Generated l incorrectly');
  end

  intersecting = 1;
  val = (u-l)/sqrt(1-u^2);
  if intersecting
    if rand < 0.5
      eps = val + rand;
    else
      eps = -val - rand;
    end
  else
    eps = -val + rand*2*val;
  end

  Q = 2*rand(n,n)-1; Q = 0.5*(Q+Q');
  c = 2*rand(n,1)-1;

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
    eigvals = sort(eig(double(Y)));
    ratio  = eigvals(n+1)/abs(eigvals(n))
%     if ratio < minratio
    if ratio < 1.0e2
%       save mysave l u eps Q c eigvals
      save(tempname('.'),'l','u','eps','Q','c','eigvals');
      minratio = ratio;
    end
%     ratios = [ratios;ratio];
  end

end

% min(ratios)
minratio
