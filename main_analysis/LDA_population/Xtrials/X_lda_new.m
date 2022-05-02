%% linear model properties
clear lmodel animalmodel
lmodel.type = 'lda';
%lmodel.boot = 0; % if 1, model returns a nonempty btstrp struct with bootstrapped vectors
lmodel.shuffle = 0; % Keep 0 to do bootstrap. If 1 or 2, instead of bootstrap it will do shuffle. 1:bsls 2:bs
lmodel.hyperparam.source = 'from_file'; %'from_params' or 'from_file': for pre-saved params
lmodel.hyperparams.delta = 0.;
lmodel.hyperparams.gamma =0.56;
lmodel.Xth =0;
lmodel.Xstyle = 2;
lmodel.nrep = 100;
lmodel.eqtrials = 1;
lmodel.doproj = 0; %0: no proj,1: autoproj(to the trials not used to fit model), 2:to a given lag
lmodel.verbose = 0;
lmodel.exptype = 'FB';
lmodel.binsizems = '150';
preSub = 0; % if 1, subtracts pre stimulus activity before doing anything else. if 2, it happens after clean remove. default is 0, 2 is for control
if preSub == 1
    preSubTxt = '_preSub';
elseif preSub == 2
    preSubTxt = '_preSub2';
else
    preSubTxt = '';
end
lmodel.filename = [lmodel.exptype,'lmodel_25Th_',lmodel.binsizems,'ms_X_style',num2str(lmodel.Xstyle),preSubTxt,'_nrep',num2str(lmodel.nrep),'.mat'];
lmodel.hyperfilename = ['hyper_',lmodel.exptype,'lmodel25Th',lmodel.binsizems,'ms_C1_R1_CV10.mat'];
lmodel.cellfiletoload = ['cells_',lmodel.exptype,'_pl20_an',lmodel.binsizems,'_lw20_exG1_onlyC1_onlyS0_plstyle1.mat'];
lmodel.appendtoanimalmodel = 1;
lmodel.minbstrialsTh = 10; % threshold for min number of baseline trials to include in analysis

%%%%%%%%%% In case params need to be overwritten from another script (used in static simulation)
% in this case, lmodel is completely overwritten and the above params are
% ineffective
if (exist('dolmodeloverwrite','var') && (dolmodeloverwrite == 1) && exist('lmodeloverwrite','var'))
    lmodel = lmodeloverwrite;
end

%% if hyper params are from file, load them
if strcmp(lmodel.hyperparam.source,'from_file') || lmodel.appendtoanimalmodel
    
    filename = fullfile('/mnt/data/Mitra/figs/P2_L/bothdircombined/2020_02_24/LDA/version2',lmodel.hyperfilename);
    h=load(filename);
end
%% cells files names and cleanup parameters (from aggregate analysis)

cleanupcriteria = struct;
cleanupcriteria.doclean = 1; % 1 or 0, perform any cleaning or not. If zero does nothing
cleanupcriteria.significance = nan; % significance of activity instead of absolute value
if  strcmp(lmodel.hyperparam.source,'from_file')
    cleanupcriteria.smb_activity = unique(cellfun(@(x) x.cleanTh, h.animalmodel)); % if not the same for all animals,
    % it errors in the next steps
else
    cleanupcriteria.smb_activity = 2.5;%
end


cellremovecriteria.doclean = 1; % 1 or 0, perform any cleaning or not. If zero does nothing
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

%% load cells, clean and remove cells: (aggregate analysis), choose animal
if ~isnan(lmodel.cellfiletoload)
    load(lmodel.cellfiletoload);
end

if preSub == 1
    [LMcells] = subtract_prestim(LMcells,lmodel.binsizems);
    [V1cells] = subtract_prestim(V1cells,lmodel.binsizems);
end

[LMcells,V1cells,Behcells,params] = cleanupcells(cleanupcriteria,LMcells,V1cells,Behcells,params);
[cell_keep_ind,LMcells,V1cells,Behcells,params] = removecells(cell_keep_ind,cellremovecriteria,LMcells,V1cells,Behcells,params);

if preSub == 2
    [LMcells] = subtract_prestim(LMcells,lmodel.binsizems);
    [V1cells] = subtract_prestim(V1cells,lmodel.binsizems);
end

% exclude animals based on minimum number of trials in baseline
numanimals = max([cellfun(@(x) x.simulcode, V1cells), cellfun(@(x) x.simulcode, LMcells)]);

includedanimals = [];
for animalnum =1:numanimals
    if strcmp(lmodel.exptype,'FF')
        targetcell = LMcells(find(cellfun(@(x) x.simulcode,LMcells) == animalnum ));
    elseif strcmp(lmodel.exptype,'FB')
        targetcell = V1cells(find(cellfun(@(x) x.simulcode,V1cells) == animalnum ));
    end
    minbstrials = unique(cellfun(@(x) min(size(x.nbs.go,1),size(x.nbs.nogo,1)),targetcell)) ;
    % should be the same number for all animals, so a scalar for minbstrials
    if minbstrials >= lmodel.minbstrialsTh
        includedanimals(end+1) = animalnum;
    end
    
end
includedanimals
numanimals = length(includedanimals);


if strcmp(lmodel.hyperparam.source,'from_file')
    nanmedian(cellfun(@(x) x.bestgamma, h.animalmodel(includedanimals)))
    nanmedian(cellfun(@(x) x.bestdelta, h.animalmodel(includedanimals)))
end

%%
nlags=8;


%%% initialization
Cmat = nan(2*nlags,2*nlags,numanimals);
Cmatweights = nan(2*nlags,2*nlags,numanimals);
Cmatdif = nan(nlags,nlags,numanimals);

godist = nan(numanimals,nlags);
nogodist = nan(numanimals,nlags);
weightgodist = nan(numanimals,nlags);
weightnogodist = nan(numanimals,nlags);

%%% combined projections from all animals
proj_bs_go = cell(1,nlags);
proj_bs_nogo = cell(1,nlags);
proj_ls_go = cell(1,nlags);
proj_ls_nogo = cell(1,nlags);
proj_dif_go = cell(1,nlags);
proj_dif_nogo = cell(1,nlags);
%
proj_animalcount_bs_go = cell(1,nlags);
proj_animalcount_bs_nogo = cell(1,nlags);
proj_animalcount_ls_go = cell(1,nlags);
proj_animalcount_ls_nogo = cell(1,nlags);

lmodel.Tefolds_go = cell(1,nlags);
lmodel.Trfolds_go = cell(1,nlags);
lmodel.Tefolds_nogo = cell(1,nlags);
lmodel.Trfolds_nogo = cell(1,nlags);

for l=1:nlags
    proj_bs_go{l} = [];
    proj_bs_nogo{l} = [];
    proj_ls_go{l} = [];
    proj_ls_nogo{l} = [];
    proj_dif_go{l} = [];
    proj_dif_nogo{l} = [];
    %
    proj_animalcount_bs_go{l} = [];
    proj_animalcount_bs_nogo{l} = [];
    proj_animalcount_ls_go{l} = [];
    proj_animalcount_ls_nogo{l} = [];
end

for repi = 1:lmodel.nrep
    lmodel.rep{repi}.part = cell(1,2);
    for part = 1:2
        lmodel.rep{repi}.part{part}.go.X = cell(1,nlags);
        lmodel.rep{repi}.part{part}.go.Y = cell(1,nlags);
        lmodel.rep{repi}.part{part}.go.Xn = cell(1,nlags);
        lmodel.rep{repi}.part{part}.go.Yn = cell(1,nlags);
        lmodel.rep{repi}.part{part}.go.vec = cell(1,nlags);
        lmodel.rep{repi}.part{part}.go.cnst = cell(1,nlags);
        lmodel.rep{repi}.part{part}.go.vec_n = cell(1,nlags);
        lmodel.rep{repi}.part{part}.go.mu_dist = cell(1,nlags);
        lmodel.rep{repi}.part{part}.go.mu_dist_norm = cell(1,nlags);
        lmodel.rep{repi}.part{part}.go.batta_dist = cell(1,nlags);
        lmodel.rep{repi}.part{part}.go.lda_norm_dist = cell(1,nlags);
        lmodel.rep{repi}.part{part}.go.alltraces = cell(1,nlags);
        
        lmodel.rep{repi}.part{part}.nogo.X = cell(1,nlags);
        lmodel.rep{repi}.part{part}.nogo.Y = cell(1,nlags);
        lmodel.rep{repi}.part{part}.nogo.Xn = cell(1,nlags);
        lmodel.rep{repi}.part{part}.nogo.Yn = cell(1,nlags);
        lmodel.rep{repi}.part{part}.nogo.vec = cell(1,nlags);
        lmodel.rep{repi}.part{part}.nogo.cnst = cell(1,nlags);
        lmodel.rep{repi}.part{part}.nogo.vec_n = cell(1,nlags);
        lmodel.rep{repi}.part{part}.nogo.mu_dist = cell(1,nlags);
        lmodel.rep{repi}.part{part}.nogo.mu_dist_norm = cell(1,nlags);
        lmodel.rep{repi}.part{part}.nogo.batta_dist = cell(1,nlags);
        lmodel.rep{repi}.part{part}.nogo.lda_norm_dist = cell(1,nlags);
        lmodel.rep{repi}.part{part}.nogo.alltraces = cell(1,nlags);
    end
end



for animalnum= includedanimals
    if strcmp(lmodel.exptype,'FF')
        targetcell = LMcells(find(cellfun(@(x) x.simulcode,LMcells) == animalnum ));
        notargetcell = V1cells(find(cellfun(@(x) x.simulcode,V1cells) == animalnum ));
    elseif strcmp(lmodel.exptype,'FB')
        targetcell = V1cells(find(cellfun(@(x) x.simulcode,V1cells) == animalnum ));
        notargetcell = LMcells(find(cellfun(@(x) x.simulcode,LMcells) == animalnum ));
    end
  
    %%% if from_file, load hyper params individually per animal
    if strcmp(lmodel.hyperparam.source,'from_file')
        lmodel.hyperparams.delta = h.animalmodel{animalnum}.bestdelta;
        lmodel.hyperparams.gamma = h.animalmodel{animalnum}.bestgamma;
    end
    %%% quick check of animal index against FBlmodel
    if strcmp(lmodel.hyperparam.source,'from_file')
        if ~(strcmp(h.animalmodel{animalnum}.animalname,unique(cellfun(@(x) x.animal, targetcell,'UniformOutput',0))))
            error('animalindices dont match')
        end
    end
    %%%
    % equalize trials
    if lmodel.eqtrials
        randeq = 0; % if zero, takes last. Seed is fixed  for random sampling
        targetcell = equalizegonogotrials(targetcell,0,0,randeq);% frnormalize
        notargetcell = equalizegonogotrials(notargetcell,0,0,randeq);% frnormalize
    end
    %%%
    % If lmodel.doproj == 2, gets a dimension to project on (lda vec of
    % proj lag). Otherwise, doesn't do anything. Just returns nan. 2
    % projection vectors: one go, and one nogo
    %     projlag=4;
    %     [proj_vec_go,proj_vec_nogo] = get_proj_vec(targetcell,lmodel,projlag);
    
    %%%
    % gonig through all lags, getting vectors, perf, etc
    for i=1:nlags
        
        for repi = 1:lmodel.nrep
            
            rng(repi,'twister')
            
            %%%% go
            % X (both parts) number of number of ncells * numtrials:
            Xbsgo = cell2mat(cellfun(@(x) x.laAbs.go{i}',targetcell,'UniformOutput',0)');
            % number of trials in Xbs, assigned 2 or 1. default in order
            Xind_bsgo = ones(size(Xbsgo,2),1);
            % Xind_bsgo(round(numel(Xind_bsgo)/2):end) = 2;
            Xind_bsgo(randsample(numel(Xind_bsgo),round(numel(Xind_bsgo)/2)))= 2;
            
            Xlsgo = cell2mat(cellfun(@(x) x.laAls.go{i}',targetcell,'UniformOutput',0)');
            Xind_lsgo = ones(size(Xlsgo,2),1);
            % Xind_lsgo(round(numel(Xind_lsgo)/2):end) = 2;
            Xind_lsgo(randsample(numel(Xind_lsgo),round(numel(Xind_lsgo)/2)))= 2;
            
            % notargetcell (only bs)\ 
            Xnbsgo = cell2mat(cellfun(@(x) x.laAbs.go{i}',notargetcell,'UniformOutput',0)');            
            Xnind_bsgo = ones(size(Xnbsgo,2),1);
            Xnind_bsgo(randsample(numel(Xnind_bsgo),round(numel(Xnind_bsgo)/2)))= 2;

            
            
            %%%% nogo
            Xbsnogo = cell2mat(cellfun(@(x) x.laAbs.nogo{i}',targetcell,'UniformOutput',0)');
            Xind_bsnogo = ones(size(Xbsnogo,2),1);
            %Xind_bsnogo(round(numel(Xind_bsnogo)/2):end) = 2;
            Xind_bsnogo(randsample(numel(Xind_bsnogo),round(numel(Xind_bsnogo)/2)))= 2;
            
            Xlsnogo = cell2mat(cellfun(@(x) x.laAls.nogo{i}',targetcell,'UniformOutput',0)');
            Xind_lsnogo = ones(size(Xlsnogo,2),1);
            %Xind_lsnogo(round(numel(Xind_lsnogo)/2):end) = 2;
            Xind_lsnogo(randsample(numel(Xind_lsnogo),round(numel(Xind_lsnogo)/2)))= 2;


            % notargetcell (only bs)\ 
            Xnbsnogo = cell2mat(cellfun(@(x) x.laAbs.nogo{i}',notargetcell,'UniformOutput',0)');            
            Xnind_bsnogo = ones(size(Xnbsnogo,2),1);
            Xnind_bsnogo(randsample(numel(Xnind_bsnogo),round(numel(Xnind_bsnogo)/2)))= 2;

            
            for part = 1:2
                if lmodel.Xstyle == 1
                    % X:(ntrialsbs+ntrialsls) * numcells,Y: (ntrialsbs+ntrialsls) * 1
                    lmodel.rep{repi}.part{part}.go.X{i} = [Xbsgo(:,find(Xind_bsgo==part))';Xlsgo(:,find(Xind_lsgo==part))'];
                    % 0: bs and 1: ls trials
                    lmodel.rep{repi}.part{part}.go.Y{i} = [zeros(numel(find(Xind_bsgo==part)),1);ones(numel(find(Xind_lsgo==part)),1)];
                    
                    bstrialind = find(Xind_bsgo==part);
                    lstrialind = find(Xind_lsgo==part);

                    % notargetcell
                    lmodel.rep{repi}.part{part}.go.Xn{i} = [Xnbsgo(:,find(Xnind_bsgo==part))'];
                    % 0: bs and 1: ls trials
                    lmodel.rep{repi}.part{part}.go.Yn{i} = [zeros(numel(find(Xnind_bsgo==part)),1)];

                elseif lmodel.Xstyle == 2
                    lmodel.rep{repi}.part{part}.go.X{i} = [Xbsgo(:,find(Xind_bsgo==part))';Xlsgo(:,:)'];
                    % 0: bs and 1: ls trials
                    lmodel.rep{repi}.part{part}.go.Y{i} = [zeros(numel(find(Xind_bsgo==part)),1);ones(size(Xlsgo(:,:),2),1)];
                    
                    bstrialind = find(Xind_bsgo==part);
                    lstrialind = 1:size(Xlsgo(:,:),2);

                    % notargetcell
                    lmodel.rep{repi}.part{part}.go.Xn{i} = [Xnbsgo(:,find(Xnind_bsgo==part))'];
                    % 0: bs and 1: ls trials
                    lmodel.rep{repi}.part{part}.go.Yn{i} = [zeros(numel(find(Xnind_bsgo==part)),1)];
                elseif lmodel.Xstyle == 0
                    lmodel.rep{repi}.part{part}.go.X{i} = [Xbsgo(:,:)';Xlsgo(:,:)'];
                    % 0: bs and 1: ls trials
                    lmodel.rep{repi}.part{part}.go.Y{i} = [zeros(size(Xbsgo(:,:),2),1);ones(size(Xlsgo(:,:),2),1)];
                    
                    bstrialind = 1:size(Xbsgo(:,:),2);
                    lstrialind = 1:size(Xlsgo(:,:),2);

                    % notargetcell
                    lmodel.rep{repi}.part{part}.go.Xn{i} = [Xnbsgo(:,:)'];
                    % 0: bs and 1: ls trials
                    lmodel.rep{repi}.part{part}.go.Yn{i} = [zeros(size(Xnbsgo(:,:),2),1)];
                    
                end
                
                % alltrials*ncells*time
                % it is the whole time, not inherently dependent on i. i
                % determines the particular selection of trials used at
                % that lag and part
                
                lmodel.rep{repi}.part{part}.go.alltraces{i} = X_make_alltraces(lmodel.rep{repi}.part{part}.go.X{i},targetcell,bstrialind,lstrialind,i,'go',lmodel);
                
                % do not manipulate X,Y befopre this point (proj depends on it)
                [lmodel.rep{repi}.part{part}.go.vec{i},lmodel.rep{repi}.part{part}.go.cnst{i},...
                    lmodel.rep{repi}.part{part}.go.vec_n{i},lmodel.rep{repi}.part{part}.go.mu_dist{i},...
                    lmodel.rep{repi}.part{part}.go.mu_dist_norm{i},lmodel.rep{repi}.part{part}.go.lda_norm_dist{i},...
                    lmodel.rep{repi}.part{part}.go.batta_dist{i}] = ...
                    X_fit_linear_classifier(lmodel,lmodel.rep{repi}.part{part}.go.X{i},lmodel.rep{repi}.part{part}.go.Y{i});
              
                
                if  lmodel.Xstyle == 1
                    lmodel.rep{repi}.part{part}.nogo.X{i} = [Xbsnogo(:,find(Xind_bsnogo==part))';Xlsnogo(:,find(Xind_lsnogo==part))'];
                    lmodel.rep{repi}.part{part}.nogo.Y{i} = [zeros(numel(find(Xind_bsnogo==part)),1);ones(numel(find(Xind_lsnogo==part)),1)];
                    
                    bstrialind = find(Xind_bsnogo==part);
                    lstrialind = find(Xind_lsnogo==part);

                    % notargetcell
                    lmodel.rep{repi}.part{part}.nogo.Xn{i} = [Xnbsnogo(:,find(Xnind_bsnogo==part))'];
                    lmodel.rep{repi}.part{part}.nogo.Yn{i} = [zeros(numel(find(Xnind_bsnogo==part)),1)];
                    
                elseif lmodel.Xstyle == 2
                    lmodel.rep{repi}.part{part}.nogo.X{i} = [Xbsnogo(:,find(Xind_bsnogo==part))';Xlsnogo(:,:)'];
                    lmodel.rep{repi}.part{part}.nogo.Y{i} = [zeros(numel(find(Xind_bsnogo==part)),1);ones(size(Xlsnogo,2),1)];
                    
                    
                    bstrialind = find(Xind_bsnogo==part);
                    lstrialind = 1:size(Xlsnogo(:,:),2);

                     % notargetcell
                    lmodel.rep{repi}.part{part}.nogo.Xn{i} = [Xnbsnogo(:,find(Xnind_bsnogo==part))'];
                    lmodel.rep{repi}.part{part}.nogo.Yn{i} = [zeros(numel(find(Xnind_bsnogo==part)),1)];
               elseif lmodel.Xstyle == 0
                    lmodel.rep{repi}.part{part}.nogo.X{i} = [Xbsnogo(:,:)';Xlsnogo(:,:)'];
                    % 0: bs and 1: ls trials
                    lmodel.rep{repi}.part{part}.nogo.Y{i} = [zeros(size(Xbsnogo(:,:),2),1);ones(size(Xlsnogo(:,:),2),1)];
                    
                    bstrialind = 1:size(Xbsnogo(:,:),2);
                    lstrialind = 1:size(Xlsnogo(:,:),2);

                    % notargetcell
                    lmodel.rep{repi}.part{part}.nogo.Xn{i} = [Xnbsnogo(:,:)'];
                    % 0: bs and 1: ls trials
                    lmodel.rep{repi}.part{part}.nogo.Yn{i} = [zeros(size(Xnbsnogo(:,:),2),1)];
                end
                
                
                lmodel.rep{repi}.part{part}.nogo.alltraces{i} = X_make_alltraces(lmodel.rep{repi}.part{part}.nogo.X{i},targetcell,bstrialind,lstrialind,i,'nogo',lmodel);
                
                [lmodel.rep{repi}.part{part}.nogo.vec{i},lmodel.rep{repi}.part{part}.nogo.cnst{i},...
                    lmodel.rep{repi}.part{part}.nogo.vec_n{i},lmodel.rep{repi}.part{part}.nogo.mu_dist{i},...
                    lmodel.rep{repi}.part{part}.nogo.mu_dist_norm{i},lmodel.rep{repi}.part{part}.nogo.lda_norm_dist{i},...
                    lmodel.rep{repi}.part{part}.nogo.batta_dist{i}] = ...
                    X_fit_linear_classifier(lmodel,lmodel.rep{repi}.part{part}.nogo.X{i},lmodel.rep{repi}.part{part}.nogo.Y{i});
               
                
            end
        end
    end
    
    if lmodel.verbose, do_plot_vecs(vec_n_go,vec_n_nogo), end
    
   
  
    %%% append to animalfile, in case you want to save afterwards
    animalmodel{animalnum}.lmodel = lmodel;

    % Cmatgo and Cmatnogo are averages over reps and parts pairwise cosine
    % similarities. dims: 8*8*nanmals (but only the relevant field full, rest nan)
    
    [Cgo,Cnogo,~,~,~,~] = X_prepareCmatrices_lda(animalmodel,animalnum);
    animalmodel{animalnum}.lmodel.Cmatgo_lda = Cgo(:,:,animalnum);
    animalmodel{animalnum}.lmodel.Cmatnogo_lda = Cnogo(:,:,animalnum);
    
    [Cgo,Cnogo,~,~,~,~] = X_prepareCmatrices_mu(animalmodel,animalnum);
    animalmodel{animalnum}.lmodel.Cmatgo_mu = Cgo(:,:,animalnum);
    animalmodel{animalnum}.lmodel.Cmatnogo_mu = Cnogo(:,:,animalnum);
    
    % target cell activity: normalized and not normalized
    [acmatgo,acmatnogo,~,~,~,~] = X_prepareCmatrices_activity(animalmodel,animalnum,lmodel,V1cells,LMcells,0);% last one: normalizaation
    animalmodel{animalnum}.lmodel.acmatgo_n0 = acmatgo(:,:,animalnum);
    animalmodel{animalnum}.lmodel.acmatnogo_n0 = acmatnogo(:,:,animalnum);
    
    [acmatgo,acmatnogo,~,~,~,~] = X_prepareCmatrices_activity(animalmodel,animalnum,lmodel,V1cells,LMcells,1);% last one: normalizaation
    animalmodel{animalnum}.lmodel.acmatgo_n1 = acmatgo(:,:,animalnum);
    animalmodel{animalnum}.lmodel.acmatnogo_n1 = acmatnogo(:,:,animalnum);
    
    % notarget cell activity: normalized and not normalized
    [acmatgo,acmatnogo,~,~,~,~] = X_prepareCmatrices_notarget_activity(animalmodel,animalnum,lmodel,V1cells,LMcells,0);% last one: normalizaation
    animalmodel{animalnum}.lmodel.acmatgo_notarget_n0 = acmatgo(:,:,animalnum);
    animalmodel{animalnum}.lmodel.acmatnogo_notarget_n0 = acmatnogo(:,:,animalnum);
    
    [acmatgo,acmatnogo,~,~,~,~] = X_prepareCmatrices_notarget_activity(animalmodel,animalnum,lmodel,V1cells,LMcells,1);% last one: normalizaation
    animalmodel{animalnum}.lmodel.acmatgo_notarget_n1 = acmatgo(:,:,animalnum);
    animalmodel{animalnum}.lmodel.acmatnogo_notarget_n1 = acmatnogo(:,:,animalnum);

   
    
end
%%%
% for each animal

[Cmatgo_lda,Cmatnogo_lda,Crepgo_lda,Crepnogo_lda,C_go_lda,C_nogo_lda] = X_prepareCmatrices_lda(animalmodel,includedanimals);
[Cmatgo_mu,Cmatnogo_mu,Crepgo_mu,Crepnogo_mu,C_go_mu,C_nogo_mu] = X_prepareCmatrices_mu(animalmodel,includedanimals);



%%%%%%%%% lda
figure;subplot(2,3,1);imagesc(nanmean(Cmatgo_lda,3));
colormap(jet(64))
set(gca,'CLim',[-0.7 0.7])
colorbar
subplot(2,3,2);imagesc(nanmean(Cmatnogo_lda,3));
colormap(jet(64))
set(gca,'CLim',[-0.7 0.7])
colorbar

s=subplot(2,3,3);
erbars_go = nan(1,8);
erbars_nogo = nan(1,8);
for i=1:8
    erbars_go(i) = 2*nanstd(C_go_lda(:,i))/sqrt(numel(find(~isnan(C_go_lda(:,i)))));
    erbars_nogo(i) = 2*nanstd(C_nogo_lda(:,i))/sqrt(numel(find(~isnan(C_nogo_lda(:,i)))));
end
errorbar([],nanmean(C_go_lda,1),erbars_go,'g')
hold on;
errorbar([],nanmean(C_nogo_lda,1),erbars_nogo,'r')
s.Title.String = 'lda';
ylim([-0.2 .6])

%%%%%%%%% mu
subplot(2,3,4);imagesc(nanmean(Cmatgo_mu,3));
colormap(jet(64))
set(gca,'CLim',[-0.7 0.7])
colorbar
subplot(2,3,5);imagesc(nanmean(Cmatnogo_mu,3));
colormap(jet(64))
set(gca,'CLim',[-0.7 0.7])
colorbar

s=subplot(2,3,6);
erbars_go = nan(1,8);
erbars_nogo = nan(1,8);
for i=1:8
    erbars_go(i) = 2*nanstd(C_go_mu(:,i))/sqrt(numel(find(~isnan(C_go_mu(:,i)))));
    erbars_nogo(i) = 2*nanstd(C_nogo_mu(:,i))/sqrt(numel(find(~isnan(C_nogo_mu(:,i)))));
end
errorbar([],nanmean(C_go_mu,1),erbars_go,'g')
hold on;
errorbar([],nanmean(C_nogo_mu,1),erbars_nogo,'r')
s.Title.String = 'mu-dif';
ylim([-.2 0.6])
%%%%%%%% consistency plots

figure;
s=subplot(2,2,1); 
histogram(C_go_lda(:,1),-0.2:0.2:1,'FaceColor','g','Normalization','pdf','EdgeAlpha',0);hold on;
histogram(C_nogo_lda(:,1),-0.2:0.2:1,'FaceColor','r','Normalization','pdf','EdgeAlpha',0);hold on;
s.Title.String = sprintf('lda-consistency-Nrep=%d,Xstyle=%d',lmodel.nrep,lmodel.Xstyle);

s=subplot(2,2,3); 
histogram(C_go_mu(:,1),-0.2:0.2:1,'FaceColor','g','Normalization','pdf','EdgeAlpha',0);hold on;
histogram(C_nogo_mu(:,1),-0.2:0.2:1,'FaceColor','r','Normalization','pdf','EdgeAlpha',0);hold on;
s.Title.String = sprintf('mu-consistency-Nrep=%d,Xstyle=%d',lmodel.nrep,lmodel.Xstyle);

s=subplot(2,2,2);
errorbar(0.9,nanmean(C_go_lda(:,1),1),2*nanstd(C_go_lda(:,1))/sqrt(numel(find(~isnan(C_go_lda(:,1))))),'g')
hold on;
errorbar(1.1,nanmean(C_nogo_lda(:,1),1),2*nanstd(C_nogo_lda(:,1))/sqrt(numel(find(~isnan(C_nogo_lda(:,1))))),'r');
hold on;
scatter(0.9,nanmean(C_go_lda(:,1),1),'.g')
hold on;
scatter(1.1,nanmean(C_nogo_lda(:,1),1),'.r');
text(0.5,0.5,sprintf('ranksum p = %d',ranksum(C_go_lda(:,1),C_nogo_lda(:,1))))
xlim([0.5 1.5])
s.Title.String = sprintf('lda-consistency-Nrep=%d,Xstyle=%d',lmodel.nrep,lmodel.Xstyle);
ylim([0 0.5])

s=subplot(2,2,4);
errorbar(0.9,nanmean(C_go_mu(:,1),1),2*nanstd(C_go_mu(:,1))/sqrt(numel(find(~isnan(C_go_mu(:,1))))),'g')
hold on;
errorbar(1.1,nanmean(C_nogo_mu(:,1),1),2*nanstd(C_nogo_mu(:,1))/sqrt(numel(find(~isnan(C_nogo_mu(:,1))))),'r');
hold on;
scatter(0.9,nanmean(C_go_mu(:,1),1),'.g')
hold on;
scatter(1.1,nanmean(C_nogo_mu(:,1),1),'.r');
xlim([0.5 1.5])
s.Title.String = sprintf('mu-consistency-Nrep=%d,Xstyle=%d',lmodel.nrep,lmodel.Xstyle);
text(0.5,0.5,sprintf('ranksum p = %d',ranksum(C_go_mu(:,1),C_nogo_mu(:,1))))
ylim([0 0.5])


%%%%%%%% saving
if lmodel.appendtoanimalmodel
    if lmodel.shuffle == 1 
        lmodel.filename = ['shuffle_',lmodel.filename];
    end
    filename = fullfile('/mnt/data/Mitra/figs/P2_L/bothdircombined/2020_02_24/LDA/version2',lmodel.filename);
    % lmodel has been added for each animal
    save(filename,'animalmodel','-v7.3');
end
