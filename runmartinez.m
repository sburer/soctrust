clear

numneed = 1000;
numgot = 0;
instance = 1;

results.cuts = [];
results.ratios = [];
results.absgaps = [];
results.relgaps1 = [];
results.relgaps2 = [];
results.secs = [];

n = 5;

while numgot < numneed

  [success,Q,c,r1,H,h,r2,cuts,eigY,LB,UB,sec] = martinez(n);

  if success
    results.cuts     = [results.cuts    ;cuts];
    results.ratios   = [results.ratios  ;eigY(n+1)/eigY(n)]; % Assumes eigY is sorted
    results.absgaps  = [results.absgaps ;UB-LB];
    results.relgaps1 = [results.relgaps1;(UB-LB)/abs(UB)];
    results.relgaps2 = [results.relgaps2;(UB-LB)/(1+abs(UB))];
    results.secs     = [results.secs    ;sec];
    if cuts > 0
      save(strcat('instance_',num2str(n),'_',num2str(instance)),'Q','c','r1','H','h','r2','cuts','eigY');
    end
    numgot = numgot+1;
    instance = instance+1;
  end

end

save(strcat('results_',num2str(n)),'results');
