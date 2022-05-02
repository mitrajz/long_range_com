% makes static effect simulatd dataset and run lda
% loads cells, cleans cells. per each basline bin, simulates equivalent
% laser trials with the mean of baseline trials in that lag, then randomly
% takes one of the lags, takes the average effect of that lag, and applies
% it to all simulated data ffrom all lags. So, if the randomly chosen
% effect is e.g. 50%, the mean effect at all 8 lags would not be exactly
% 50, but fluctuating around it, reflecting the noise in effcet that would
% have been cause by a constant underlying effect over time + poisson
% sampling noise.

% after simulated data has been made, it is passed to lda for com vectors
% to be calculated. The results are in cntrl, which need to be saved
% manually.
% Last section takes a very long time for plotting (although texts produced earlier),
% terminate before ending if no hist plots neede. 
% run this, save cntrl, then run X_tc_staticeffectsimulation.m

%%
% params:


nran = 100;% 100
exptype = 'FB';
nlda = 100; % 100: no other lda param needs changing

% Bellow crits are for static hist criteria and simulation sampling. Once the simulated data is fed to lda,
% it is processed with the lda criteria (clean1 and remove1)
cleancrit = 1;
removecrit = 0;

calcCom = 1; % keep 1 if you want to actually run the lda code and calculate communication.
% 0 is only used to make static hists only (calclatin allV1 and allLM)

% save mannually (commented)

%%
clear lmodel animalmodel
lmodel.type = 'lda';

lmodel.shuffle = 0; % Keep 0 to do bootstrap. If 1 or 2, instead of bootstrap it will do shuffle. 1:bsls 2:bs
lmodel.hyperparam.source = 'from_file'; %'from_params' or 'from_file': for pre-saved params
lmodel.hyperparams.delta = 0.;
lmodel.hyperparams.gamma =0.56;
lmodel.Xth = 0; % def : 0 

lmodel.Xstyle = 2;% 
lmodel.nrep = nlda;% 100
lmodel.eqtrials = 1;
lmodel.doproj = 0; %0: no proj,1: autoproj(to the trials not used to fit model), 2:to a given lag
lmodel.verbose = 0;
lmodel.exptype = exptype;
lmodel.binsizems = '150';
lmodel.filename = [lmodel.exptype,'lmodel_25Th_',lmodel.binsizems,'ms_X_style',num2str(lmodel.Xstyle),'_nrep',num2str(lmodel.nrep),'.mat'];
lmodel.hyperfilename = ['hyper_',lmodel.exptype,'lmodel25Th',lmodel.binsizems,'ms_C1_R1_CV10.mat'];
lmodel.cellfiletoload = nan;
lmodel.appendtoanimalmodel = 0;
lmodel.minbstrialsTh = 10;

%%%%%%%%%%%% inorder to overwrite X_lda_new params
lmodeloverwrite = lmodel;
dolmodeloverwrite = 1;
%% make static plots for the original dataset
load(['cells_',exptype,'_pl20_an150_lw20_exG1_onlyC1_onlyS0_plstyle1.mat']);
%%%%%%
cleanupcriteria = struct;
cleanupcriteria.doclean = cleancrit; % 1 or 0, perform any cleaning or not. If zero does nothing
cleanupcriteria.significance = nan; % significance of activity instead of absolute value
cleanupcriteria.smb_activity = 2.5; % in Hz, 

cellremovecriteria.doclean = removecrit; % 1 or 0, perform any cleaning or not. If zero does nothing
cell_keep_ind = nan; % important param: if given, indexes are no recalculated
cellremovecriteria.stability = 0;
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

%%% clean and remove cells

cell_keep_ind = nan; % important param: if given, indexes are no recalculated
[LMcells,V1cells,Behcells,params] = cleanupcells(cleanupcriteria,LMcells,V1cells,Behcells,params);
[cell_keep_ind,LMcells,V1cells,Behcells,params] = removecells(cell_keep_ind,cellremovecriteria,LMcells,V1cells,Behcells,params);

%%%%%

globalplot.histnorm = 'pdf';
if strcmp(exptype,'FB')
    targetcellname = 'V1cells';
elseif strcmp(exptype,'FF')
    targetcellname = 'LMcells';
end
subselectcells = 0;
static_effectsize_histograms;

%% assumption is that %change is constant not delta:
% IMPORTANT: set X_lda_new.m initial parameters
allV1_sim = [];
allLM_sim = [];
cntrl = cell(1,nran);

Behcells = [];
method = 1;
rng_used = 1;


for ran = 1:nran
    ran
    
    load(['cells_',exptype,'_pl20_an150_lw20_exG1_onlyC1_onlyS0_plstyle1.mat']);
    
    cleanupcriteria.doclean = cleancrit; %1
    cellremovecriteria.doclean = removecrit;% 0
    cell_keep_ind = nan; % important param: if given, indexes are no recalculated
    [LMcells,V1cells,Behcells,params] = cleanupcells(cleanupcriteria,LMcells,V1cells,Behcells,params);
    [cell_keep_ind,LMcells,V1cells,Behcells,params] = removecells(cell_keep_ind,cellremovecriteria,LMcells,V1cells,Behcells,params);
    
    if strcmp(exptype,'FB')
        targetcell = V1cells;
    elseif strcmp(exptype,'FF')
        targetcell = LMcells;
    end
    
    rng(ran,'twister');
    lagorder = randi(8,2,numel(targetcell));
    
    
    % only changing target cells
    for cellid = 1: numel(targetcell)
        
        randlag = lagorder(:,cellid);
        effect_go = targetcell{cellid}.smb.go(randlag(1));
        effect_nogo = targetcell{cellid}.smb.nogo(randlag(2));
        
        if method == 2
            smbdistgo = cell2mat( cellfun(@(x,y) 100*(x - nanmean(y))/nanmean(y),...
                targetcell{cellid}.laAls.go,targetcell{cellid}.laAbs.go,'UniformOutput',0)' );
            smbdistgo(isinf(smbdistgo)) = nan;
            try % 
                pd_go = fitdist(smbdistgo,'Kernel','Width',1);
            catch
                pd_go = [];
            end
            
            smbdistnogo = cell2mat( cellfun(@(x,y) 100*(x - nanmean(y))/nanmean(y),...
                targetcell{cellid}.laAls.nogo,targetcell{cellid}.laAbs.nogo,'UniformOutput',0)' );
            
            smbdistnogo(isinf(smbdistnogo)) = nan;
            try
                pd_nogo = fitdist(smbdistnogo,'Kernel','Width',1);
            catch
                pd_nogo = [];
            end
        end
        
        for lag = 1:8
            
            %%%% same number of samples for go and nogo
            lag_nsamples = min(numel(targetcell{cellid}.laAls.go{lag}),numel(targetcell{cellid}.laAls.nogo{lag}));
            %%%%% go
            %%%%% assign laAls
            if method == 1
                % trial averaged activity of ls trials, calculated from replaced
                % effect
                mean_bs = nanmean(targetcell{cellid}.laAbs.go{lag});
                % make replacement activity from poisson dist, with mean_ls as lambda
                rng(rng_used,'twister');
                rng_used = rng_used+1;
                replacement = poissrnd(mean_bs,numel(targetcell{cellid}.laAls.go{lag}),1);
                targetcell{cellid}.laAls.go{lag} = (1+effect_go/100) * replacement;
            elseif method == 2
                if isempty(pd_go)
                    targetcell{cellid}.laAls.go{lag} = nan(size(targetcell{cellid}.laAls.go{lag}));
                else
                    rng(rng_used,'twister');
                    rng_used = rng_used+1;
                    effects = random(pd_go,1,numel(targetcell{cellid}.laAls.go{lag}))';
                    targetcell{cellid}.laAls.go{lag} = (1+effects/100) * nanmean(targetcell{cellid}.laAbs.go{lag});
                end
            end
            
            targetcell{cellid}.smb.go(lag) = 100 * (nanmean(targetcell{cellid}.laAls.go{lag}) - nanmean(targetcell{cellid}.laAbs.go{lag}))/...
                nanmean(targetcell{cellid}.laAbs.go{lag});
            
            
            %%%%% nogo
            if method == 1
                mean_bs = nanmean(targetcell{cellid}.laAbs.nogo{lag});
                % make replacement activity from poisson dist, with mean_ls as lambda
                rng(rng_used,'twister');
                rng_used = rng_used+1;
                replacement = poissrnd(mean_bs,numel(targetcell{cellid}.laAls.nogo{lag}),1);
                targetcell{cellid}.laAls.nogo{lag} = (1+effect_nogo/100) * replacement;
            elseif method == 2
                if isempty(pd_nogo)
                    targetcell{cellid}.laAls.nogo{lag} = nan(size(targetcell{cellid}.laAls.nogo{lag}));
                else
                    rng(rng_used,'twister');
                    rng_used = rng_used+1;
                    effects = random(pd_nogo,1,numel(targetcell{cellid}.laAls.nogo{lag}))';
                    targetcell{cellid}.laAls.nogo{lag} = (1+effects/100) * nanmean(targetcell{cellid}.laAbs.nogo{lag});
                end
            end
            
            targetcell{cellid}.smb.nogo(lag) = 100 * (nanmean(targetcell{cellid}.laAls.nogo{lag}) - nanmean(targetcell{cellid}.laAbs.nogo{lag}))/...
                nanmean(targetcell{cellid}.laAbs.nogo{lag});
        end
        
    end
    
    if strcmp(exptype,'FB')
        V1cells = targetcell;
    elseif strcmp(exptype,'FF')
        LMcells = targetcell;
    end
    
    allV1_sim = [allV1_sim,V1cells];
    allLM_sim = [allLM_sim,LMcells];
    
    
    if calcCom
        run('X_lda_new.m')
        cntrl{ran}.go = cellfun(@(x) x.lmodel.Cmatgo_lda,animalmodel,'UniformOutput',0);
        cntrl{ran}.nogo = cellfun(@(x) x.lmodel.Cmatnogo_lda,animalmodel,'UniformOutput',0);
        close all
        clear animalmodel
    end
    
end


%save('cntrl_ff_new_method1_nomeanadj_nran100_lda100_R0C1.mat','cntrl','-v7.3')

%%
V1cells =  allV1_sim;
LMcells =  allLM_sim;

globalplot.histnorm = 'pdf';
if strcmp(exptype,'FB')
    targetcellname = 'V1cells';
elseif strcmp(exptype,'FF')
    targetcellname = 'LMcells';
end
subselectcells = 0;
static_effectsize_histograms;

