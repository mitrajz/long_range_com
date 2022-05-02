% saves lag0,lg1 and lag8 per animal and lag. no fit params saved, fits
% dont matter at all here
% Do the following:
% open X_population_pc_temporal_drift.m
% comment out included animals: 
% comment out the entire lsbs for loop

% gopcang has 8 values that are lags 0 - 7
% there is never a lag8! because lag is delta t and is maximum 7!

savdir = '/mnt/data/Mitra/figs/P2_L/bothdircombined/2020_02_24/LDA/version2/cosinesimilaritymats_and_timeconstants/saveres';
saveres = 1;


animal = cell(1,6);
for animali= 1:6

    animal{animali}.pcn = cell(1,3);
    includedanimals = animali;
    
    
    lsbs = 1;
    X_population_pc_temporal_drift;
    
    for i=1:3
        animal{animali}.pcn{i}.ls.go.lag0 = gopcang{1}(2:end,i);
        animal{animali}.pcn{i}.ls.go.lag1 = [gopcang{2}(2:end,i);nan];
        animal{animali}.pcn{i}.ls.go.lag8 = [gopcang{8}(2:end,i);nan(7,1)];

        animal{animali}.pcn{i}.ls.nogo.lag0 = nogopcang{1}(2:end,i);
        animal{animali}.pcn{i}.ls.nogo.lag1 = [nogopcang{2}(2:end,i);nan];
        animal{animali}.pcn{i}.ls.nogo.lag8 = [nogopcang{8}(2:end,i);nan(7,1)];
    end
    
    lsbs = 2;
    X_population_pc_temporal_drift;

    for i=1:3
        animal{animali}.pcn{i}.bs.go.lag0 = gopcang{1}(2:end,i);
        animal{animali}.pcn{i}.bs.go.lag1 = [gopcang{2}(2:end,i);nan];
        animal{animali}.pcn{i}.bs.go.lag8 = [gopcang{8}(2:end,i);nan(7,1)];

        animal{animali}.pcn{i}.bs.nogo.lag0 = nogopcang{1}(2:end,i);
        animal{animali}.pcn{i}.bs.nogo.lag1 = [nogopcang{2}(2:end,i);nan];
        animal{animali}.pcn{i}.bs.nogo.lag8 = [nogopcang{8}(2:end,i);nan(7,1)];
    end
   

end

if saveres
    cd(savdir)
    foldername = ['pc_',exptype,binsizems,'Xstyle',xstyle,'nrep',num2str(nrep),'_','eqbsls',num2str(eqbsls),...
        'numtrials',num2str(numtrials),'_',datestr(now,'yyyy_mm_dd_HH_MM_SS')];
    mkdir(foldername)
    cd(foldername)
    % figures
   
    % variables
    save('allanimalpcs.mat','animal')
end
