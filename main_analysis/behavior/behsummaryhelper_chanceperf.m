function    [dp] = behsummaryhelper_chanceperf(gotrialind,nogotrialind,...
        correctgotrialind,incorrectgotrialind,correctnogotrialind,incorrectnogotrialind,...
        nogroomingind,onlynogrooming);

nrep  = 5000;

if onlynogrooming
    correctgotrialind = intersect(correctgotrialind,nogroomingind);
    incorrectgotrialind = intersect(incorrectgotrialind,nogroomingind);
    correctnogotrialind = intersect(correctnogotrialind,nogroomingind);
    incorrectnogotrialind = intersect(incorrectnogotrialind,nogroomingind);
    
    gotrialind = intersect(gotrialind,nogroomingind);
    nogotrialind = intersect(nogotrialind,nogroomingind);
end

licked = [correctgotrialind incorrectnogotrialind];
nolicked = [correctnogotrialind incorrectgotrialind];


pooledtrials = [gotrialind' nogotrialind'];

dp = nan(1,nrep);

for rep = 1:nrep
    pooledtrialoutcome = nan(size(pooledtrials));
    
    pooledrandlick = randi(size(pooledtrials),[1,length(licked)+length(nolicked)]);
    pooledtrialoutcome(pooledrandlick(1:length(licked))) = 1;
    pooledtrialoutcome(pooledrandlick((length(licked)+1):end)) = 0;
    
    perm_correctgotrialind = gotrialind( find(pooledtrialoutcome(1:length(gotrialind)) == 1) );
    perm_incorrectgotrialind = gotrialind( find(pooledtrialoutcome(1:length(gotrialind)) == 0) );
    perm_correctnogotrialind = nogotrialind( find(pooledtrialoutcome((1+length(gotrialind)):end) == 0) );
    perm_incorrectnogotrialind = nogotrialind( find(pooledtrialoutcome((1+length(gotrialind)):end) == 1) );
    
    % this doesn't count early lick trials
    H = length(perm_correctgotrialind)/(length(perm_correctgotrialind)+length(perm_incorrectgotrialind));
    
    F = length(perm_incorrectnogotrialind)/(length(perm_correctnogotrialind)+length(perm_incorrectnogotrialind));
    
    dp(rep) = norminv(H) - norminv(F);    
end