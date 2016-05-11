n = 3;

ratios = [];

minratio = +Inf;

for iter = 1:1


     Q = [ 0.462352790165023   0.119315430204378   0.634836483510965    % Has about a 1% gap
           0.119315430204378  -0.942540348543900  -0.408067614985745
           0.634836483510965  -0.408067614985745   0.086272852556693 ]; 

           c = [-0.250429143407273
                 0.887453125996381
                -0.773145049597823];

           l = -0.201680972975058;
           u =  0.267392086634794;

           eps = 0.856586099013129;


     Q = [ 0.46   0.11   0.63    % Has about a 1% gap
           0.11  -0.94  -0.40
           0.63  -0.40   0.08 ];

           c = [-0.25
                 0.88
                -0.77];

           l = -0.20;
           u =  0.26;

           eps = 0.85;

%      Q = [ 0.4   0.1   0.6
%            0.1  -0.9  -0.4
%            0.6  -0.4   0.0 ];

%            c = [-0.2
%                  0.8
%                 -0.7];

%            l = -0.2;
%            u =  0.2;

%            eps = 0.8;

%      Q = [ 4   0   6
%            0  -8  -4
%            6  -4   0 ];

%            c = [-2
%                  8
%                 -6];

%            l = -0.2;
%            u =  0.2;

%            eps = 0.7;


%      Q = [ 2   0   3 
%            0  -4  -2
%            3  -2   0 ];

%            c = [-1
%                  4 
%                 -3];

%            l = -0.2;
%            u =  0.2;

%            eps = 0.7;


%      Q = [ 3   1   3     % This one has about 0.6% gap
%            1  -4  -2
%            3  -2   0];

%            c = [-1
%                  4
%                 -3];

%            l = -0.2;
%            u =  0.2;

%            eps = 0.7;


% With eps = 0.7, two cuts appear to intersect inside circle.

%   l = -rand;
%   u =  rand;

%   eps = 2*rand(1)-1;

%   Q = 2*rand(n,n)-1; Q = 0.5*(Q+Q');
%   c = 2*rand(n,1)-1;

  load('mytmp.mat');
%   load('tp232270299193482.mat');

% Q = [ 0.119602205174211   0.148272624702926   0.589211451155130
%       0.148272624702926  -0.932905803142221   0.316686456557858
%       0.589211451155130   0.316686456557858   0.012850573747279 ];

% c = [ 0.698614915960597
%       0.747151491137563
%       0.440821720650267 ];

% l = -0.433654760497416;

% u =  0.027915208058777;

% eps = 1.242036271326532;

% Q = [ 0.1196   0.1482   0.5892
%       0.1482  -0.9329   0.3166
%       0.5892   0.3166   0.0128 ];

% c = [ 0.6986
%       0.7471
%       0.4408 ];

% l = -0.4336;

% u =  0.0279;

% eps = 1.2420;

% Q = [ 0.11   0.14   0.58
%       0.14  -0.93   0.31
%       0.58   0.31   0.01 ];

% c = [ 0.69
%       0.74
%       0.44 ];

% l = -0.43;

% u =  0.02;

% eps = 1.24;

Q = [ 0.10   0.15   0.60
      0.15  -0.95   0.30
      0.60   0.30   0.00 ];

c = [ 0.70
      0.70
      0.45 ];

Q = 100*Q/5
c = 100*c/5

l = -0.5

u =  0.0

eps = 1.2

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
    fprintf('keep %.12f\n', double(obj));
    eigvals = sort(eig(double(Y)))
    ratio  = eigvals(n+1)/abs(eigvals(n));
    if ratio < minratio
%       save mysave l u eps Q c eigvals
      minratio = ratio;
    end
    ratios = [ratios;ratio];
  end

end

min(ratios)


fid = fopen('trs2_fixed.dat','w');

fprintf(fid,'param Q:    1   2   3 :=\n');
for i = 1:n
  fprintf(fid,'%d  ', i);
  for j = 1:n
    fprintf(fid,'%.15f  ', Q(i,j));
  end
  fprintf(fid,'\n');
end
fprintf(fid,';\n');

fprintf(fid,'param c :=\n');
for j = 1:n
  fprintf(fid,'%d   %.15f\n', j, c(j));
end
fprintf(fid,';\n');

fprintf(fid,'param l := %.15f;\n', l);
fprintf(fid,'param u := %.15f;\n', u);
fprintf(fid,'param eps := %.15f;\n', eps);

fclose(fid);

% unix('ampl < trs2_fixed.ampl');

