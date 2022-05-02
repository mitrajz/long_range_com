mean_standarderr_coeff =2.73%1.96; % % Bonferroni corrected 2 sided comparison: (1-normcdf(1.96))/(1-normcdf(2.73)) = 8
nrep = 100;% 100
%% first load cells , only clean (according to latest methods) R0C1
exptype = 'FB';
load(['cells_',exptype,'_pl20_an150_lw20_exG1_onlyC1_onlyS0_plstyle1.mat']);


cleanupcriteria = struct;
cleanupcriteria.doclean = 1; % 1 or 0, perform any cleaning or not. If zero does nothing
cleanupcriteria.significance = nan; % significance of activity instead of absolute value
cleanupcriteria.smb_activity = 2.5; % in Hz

cellremovecriteria.doclean = 0; % 1 or 0, perform any cleaning or not. If zero does nothing
cell_keep_ind = nan; % important param: if given, indexes are no recalculated
cellremovecriteria.stability = 1;
cellremovecriteria.lindriftTH = 0.5; % stability param - def:0.5
cellremovecriteria.activity = nan;
cellremovecriteria.responsiveness =[4 1];  % nan or several options:
cellremovecriteria.sigTH = 0.05;
cellremovecriteria.minspkTH = 1; %responsiveness param (method 4 only)- def:1
cellremovecriteria.showexampleplots = 0; % order:targetgo, targetnogo, not_target go, not_targetnogo

[LMcells,V1cells,Behcells,params] = cleanupcells(cleanupcriteria,LMcells,V1cells,Behcells,params);
[cell_keep_ind,LMcells,V1cells,Behcells,params] = removecells(cell_keep_ind,cellremovecriteria,LMcells,V1cells,Behcells,params);


%%
howmanysigbins = nan(max( unique(cellfun(@(x) x.simulcode, LMcells))),9); % here Lm or V1 doent matter
howmanycellssig = nan(max( unique(cellfun(@(x) x.simulcode, LMcells))),8);
for animaln = unique(cellfun(@(x) x.simulcode, LMcells))
    if strcmp(exptype,'FF')
        targetcell = LMcells(find(cellfun(@(x) x.simulcode, LMcells) == animaln));
    elseif strcmp(exptype,'FB')
        targetcell = V1cells(find(cellfun(@(x) x.simulcode, V1cells) == animaln));
    end
    
    allgo = cell2mat(cellfun(@(x) x.smb.go, targetcell,'UniformOutput',0));
    allnogo = cell2mat(cellfun(@(x) x.smb.nogo, targetcell,'UniformOutput',0));
    all = cell2mat([cellfun(@(x) x.smb.nogo, targetcell,'UniformOutput',0),cellfun(@(x) x.smb.go, targetcell,'UniformOutput',0)]);
    
    % standard error of the mean
    standarderr_prct_go = nan(length(targetcell),length(targetcell{1}.laAbs.go));
    standarderr_prct_nogo = nan(length(targetcell),length(targetcell{1}.laAbs.nogo));
    for i=1:length(targetcell)
        % go
        for l=1:length(targetcell{i}.laAbs.go)
            mean_standarderr = nanstd(targetcell{i}.laAbs.go{l})/sqrt(numel(targetcell{i}.laAbs.go{l}));
            standarderr_prct_go(i,l) = 100 * mean_standarderr_coeff * mean_standarderr/nanmean(targetcell{i}.laAbs.go{l});
            if isinf(standarderr_prct_go(i,l))
                standarderr_prct_go(i,l) = nan;
            end
        end
        % nogo
        for l=1:length(targetcell{i}.laAbs.nogo)
            mean_standarderr = nanstd(targetcell{i}.laAbs.nogo{l})/sqrt(numel(targetcell{i}.laAbs.nogo{l}));
            standarderr_prct_nogo(i,l) = 100 * mean_standarderr_coeff * mean_standarderr/nanmean(targetcell{i}.laAbs.nogo{l});
            if isinf(standarderr_prct_nogo(i,l))
                standarderr_prct_nogo(i,l) = nan;
            end
        end
    end
    
   
    allmeanconfs = [reshape(standarderr_prct_nogo',1,[]),reshape(standarderr_prct_go',1,[])];
    gomeanconfs = reshape(standarderr_prct_go',1,[]);
    nogomeanconfs = reshape(standarderr_prct_nogo',1,[]);
    
    all_s = all; allgo_s = allgo; allnogo_s = allnogo;
    all_s(abs(all_s)<allmeanconfs) = nan;
    allgo_s(abs(allgo_s)<gomeanconfs) = nan;
    allnogo_s(abs(allnogo_s)<nogomeanconfs) = nan;
    
    % this is how many cells were affected at least at one point.
    % color correct the pie charts across ff and fb posthoc
    % reshape(allgo_s,8,[]) is 8*numcells. nan for non sig effects, and the
    % effsct percentage (smb) for rest
    
    nsigbins_go = 100*[numel(find(sum(isnan(reshape(allgo_s,8,[])),1)==8))/numel(targetcell),...
        numel(find(sum(isnan(reshape(allgo_s,8,[])),1)==7))/numel(targetcell),...
        numel(find(sum(isnan(reshape(allgo_s,8,[])),1)==6))/numel(targetcell),...
        numel(find(sum(isnan(reshape(allgo_s,8,[])),1)==5))/numel(targetcell),...
        numel(find(sum(isnan(reshape(allgo_s,8,[])),1)==4))/numel(targetcell),...
        numel(find(sum(isnan(reshape(allgo_s,8,[])),1)==3))/numel(targetcell),...
        numel(find(sum(isnan(reshape(allgo_s,8,[])),1)==2))/numel(targetcell),...
        numel(find(sum(isnan(reshape(allgo_s,8,[])),1)==1))/numel(targetcell),...
        numel(find(sum(isnan(reshape(allgo_s,8,[])),1)==0))/numel(targetcell)];
    atleast1_go = 100*numel(find(sum(isnan(reshape(allgo_s,8,[])),1)<8))/numel(targetcell);
    
    nsigbins_nogo = 100*[numel(find(sum(isnan(reshape(allnogo_s,8,[])),1)==8))/numel(targetcell),...
        numel(find(sum(isnan(reshape(allnogo_s,8,[])),1)==7))/numel(targetcell),...
        numel(find(sum(isnan(reshape(allnogo_s,8,[])),1)==6))/numel(targetcell),...
        numel(find(sum(isnan(reshape(allnogo_s,8,[])),1)==5))/numel(targetcell),...
        numel(find(sum(isnan(reshape(allnogo_s,8,[])),1)==4))/numel(targetcell),...
        numel(find(sum(isnan(reshape(allnogo_s,8,[])),1)==3))/numel(targetcell),...
        numel(find(sum(isnan(reshape(allnogo_s,8,[])),1)==2))/numel(targetcell),...
        numel(find(sum(isnan(reshape(allnogo_s,8,[])),1)==1))/numel(targetcell),...
        numel(find(sum(isnan(reshape(allnogo_s,8,[])),1)==0))/numel(targetcell)];
    atleast1_nogo = 100*numel(find(sum(isnan(reshape(allnogo_s,8,[])),1)<8))/numel(targetcell);
    
    nsigbins_av = nanmean([nsigbins_go;nsigbins_nogo]); % this is average over go and nogo for the1 animal
    howmanysigbins(animaln,:) = nsigbins_av;
    %atleast1_av = nanmean([atleast1_go,atleast1_nogo])
    
    %%%%%%%%% How many cells sig per bin
    %This is an 8*2 matrix. 8 bins(percentage of cells sig in that bin) and the
    %2 column are go and nogo. Now these 2 column could be averaged or
    %sumed
    bothgng = 100*[sum(~isnan(reshape(allgo_s,8,[])),2)./(size(allgo_s,2)/8)...
        ,sum(~isnan(reshape(allnogo_s,8,[])),2)./(size(allnogo_s,2)/8)];
    
    howmanycellssig(animaln,:) = nanmean(bothgng,2);
end
figure
for l=1:9
    if l == 1
        bar(l-1,nanmean(howmanysigbins(:,l)),'EdgeAlpha',0,'FaceColor',[0.5 0.5 0.5])
    else
        bar(l-1,nanmean(howmanysigbins(:,l)),'EdgeAlpha',0,'FaceColor','b')
    end
    hold on;errorbar(l-1,nanmean(howmanysigbins(:,l)),...
    2*nanstd(howmanysigbins(:,l))/sqrt(numel(howmanysigbins(:,l))),...
    'color','k')
end
title(exptype);ylabel('percentage of cells');xlabel('number of time windows with significant effect');
set(gca,'XTick',0:8)
%%% the number of cells affected in at least 1 timebin
[exptype,' - ' , sprintf('percentage of cells sig affected in at least 1 time bin = %f +- %f (mean+-s.e.m across animals)',...
    nanmean(100-howmanysigbins(:,1)),2*nanstd(100-howmanysigbins(:,1))./sqrt(numel(howmanycellssig(:,1))) )]

figure
for l=1:8
    bar(l,nanmean(howmanycellssig(:,l)),'EdgeAlpha',0,'FaceColor','b')

    hold on;errorbar(l,nanmean(howmanycellssig(:,l)),...
    2*nanstd(howmanycellssig(:,l))/sqrt(numel(howmanycellssig(:,l))),...
    'color','k')
end
title([exptype,' - f-test(anova) p = ',num2str(anova1(howmanycellssig,[],'off'))]);ylabel('percentage of cells');xlabel('silencing time window');
set(gca,'XTick',1:8);set(gca,'XTickLabel',{'T1','T2','T3','T4','T5','T6','T7','T8'})
