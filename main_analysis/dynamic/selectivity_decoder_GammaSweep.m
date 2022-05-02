
% This is to be run before selectivity_decoder.m to get best gammas.
% importantly nfold, nrep and the trial selection
% should be identical in the two. 

Exptype = 'FF';


cleanupcriteria = struct;
cleanupcriteria.doclean = 1; % 1 or 0, perform any cleaning or not. If zero does nothing
cleanupcriteria.significance = nan; % significance of activity instead of absolute value
cleanupcriteria.smb_activity = 2.5; % in Hz,

cellremovecriteria.doclean = 0; % 1 or 0, perform any cleaning or not. If zero does nothing
cell_keep_ind = nan; % important param: if given, indexes are no recalculated
cellremovecriteria.stability = 1;
cellremovecriteria.lindriftTH = 0.5; % stability param - def:0.5
cellremovecriteria.activity = nan;
cellremovecriteria.responsiveness =[4 1];  % nan or several options:
% first element: method
% second element: at least one of
% go or nogo, or both: 1 or 2
% def: 4,1
cellremovecriteria.sigTH = 0.05; % responsiveness param - def:0.05
cellremovecriteria.minspkTH = 1; %responsiveness param (method 4 only)- def:1
cellremovecriteria.showexampleplots = 0; % order:targetgo, targetnogo, not_target go, not_targetnogo


if strcmp(Exptype,'FB')
    load('cells_FB_pl20_an80_lw20_exG1_onlyC1_onlyS0_plstyle1.mat');
elseif strcmp(Exptype,'FF')
    load('cells_FF_pl20_an80_lw20_exG1_onlyC1_onlyS0_plstyle1.mat');
end

[LMcells,V1cells,Behcells,params] = cleanupcells(cleanupcriteria,LMcells,V1cells,Behcells,params);
[cell_keep_ind,LMcells,V1cells,Behcells,params] = removecells(cell_keep_ind,cellremovecriteria,LMcells,V1cells,Behcells,params);


%%% CV with folds, and add time axis (steps of mean and ste). Add pvalues

clc


nfold = 5;
nrep = 20;
if strcmp(Exptype,'FB')
    numanimals = 6;
elseif strcmp(Exptype,'FF')
    numanimals = 7;
end
timeWindowSec = 0.08;
gammarange = 0:.1:1;
ngamma = numel(gammarange);

numtimepoints = 8;


res_mat = nan(numanimals,ngamma);


for gamma = 1:ngamma
    bs_mat = nan(numanimals,numtimepoints);
    for timepoint = 1:numtimepoints
        
        ClAc_bs = nan(numanimals,nrep);
        
        for animalnum=1:numanimals
            for rep=1:nrep
                
                if strcmp(Exptype,'FB')
                    targetcell = V1cells(find(cellfun(@(x) x.simulcode,V1cells) == animalnum ));
                elseif strcmp(Exptype,'FF')
                    targetcell = LMcells(find(cellfun(@(x) x.simulcode,LMcells) == animalnum ));
                end
                % ntrials*ncells:
                gobs = cell2mat(cellfun(@(x) x.laAbs.go{timepoint},targetcell,'UniformOutput',0));
                nogobs = cell2mat(cellfun(@(x) x.laAbs.nogo{timepoint},targetcell,'UniformOutput',0));
                gols = cell2mat(cellfun(@(x) x.laAls.go{timepoint},targetcell,'UniformOutput',0));
                nogols = cell2mat(cellfun(@(x) x.laAls.nogo{timepoint},targetcell,'UniformOutput',0));
                
                bsall = [gobs;nogobs];
                lsall = [gols;nogols];
                
                
                % equalize trials
                gobs(:,isnan(sum(bsall,1))) = [];
                nogobs(:,isnan(sum(bsall,1))) = [];
                gols(:,isnan(sum(lsall,1))) = [];
                nogols(:,isnan(sum(lsall,1))) = [];
                
                mintrialnum = min([size(gobs,1),size(nogobs,1),size(gols,1),size(nogols,1)]);
                gobs = gobs(1:mintrialnum,:);
                nogobs = nogobs(1:mintrialnum,:);
                
                %
                rng(rep)

                try

                    cvind=crossvalind('Kfold',size(gobs,1)+size(nogobs,1),nfold);
                    ClAc_bs_fold = [];
                    
                    for fold=1:nfold
                        traininds = (find(cvind ~= fold));
                        testinds = (find(cvind == fold));
                        
                        PooledX = [gobs;nogobs];
                        PooledY = [zeros(size(gobs,1),1);ones(size(nogobs,1),1)];
                        gammarange(gamma)
                        mdl_bs = fitcdiscr(PooledX(traininds,:),...
                            PooledY(traininds,:),'DiscrimType','linear','Prior','uniform','Gamma',gammarange(gamma));
                        predicted_label = predict(mdl_bs,PooledX(testinds,:));
                        real_label = PooledY(testinds);
                        % check classification accuracy
                        ClAc_bs_fold(end+1) = sum(~xor(predicted_label,real_label))/numel(predicted_label);%nanmean((nanmean(gobs,1)+nanmean(nogobs,1))/2);%
                        
                        
                    end
                    ClAc_bs(animalnum,rep) = nanmean(ClAc_bs_fold);
                    
                catch
                end
            end
            bs_mat(animalnum,timepoint) = nanmean(ClAc_bs(animalnum,:));
            
        end
    end
    res_mat(:,gamma) = nanmean(bs_mat,2);
end
%%
[m,i]=max(res_mat')
gammarange(i)