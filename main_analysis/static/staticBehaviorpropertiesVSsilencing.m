% Before running this script, make sure retin_pipeline.m has been run for
% the animals, depth and spikeshape has been assiged, and
% makecells_generic* has been run for all included animals with retin,
% depth, and spikeshape assignment
%% loading and preparing

exptype = 'FF';

cleanupcriteria = struct;
cleanupcriteria.doclean = 1; % 1 or 0, perform any cleaning or not. If zero does nothing
cleanupcriteria.significance = nan; % significance of activity instead of absolute value
cleanupcriteria.smb_activity = 2.5; % in Hz
cellremovecriteria.doclean = 0; % 1 or 0, perform any cleaning or not. If zero does nothing
cell_keep_ind = nan; % important param: if given, indexes are no recalculated
cellremovecriteria.stability = 1;
cellremovecriteria.lindriftTH = 0.5; % stability param - def:0.5
cellremovecriteria.activity = nan;
cellremovecriteria.responsiveness =[4 1];
cellremovecriteria.sigTH = 0.05; % responsiveness param - def:0.05
cellremovecriteria.minspkTH = 1; %responsiveness param (method 4 only)- def:1
cellremovecriteria.showexampleplots = 0; % order:targetgo, targetnogo, not_target go, not_targetnogo

if strcmp(exptype,'FF')
    C = load('cells_FF_pl20_an150_lw20_exG1_onlyC1_onlyS0_plstyle1.mat');
    IC = load('cells_FF_pl20_an150_lw20_exG0_onlyC-1_onlyS0_plstyle1.mat'); 
else
    C = load('cells_FB_pl20_an150_lw20_exG1_onlyC1_onlyS0_plstyle1.mat');
    IC = load('cells_FB_pl20_an150_lw20_exG0_onlyC-1_onlyS0_plstyle1.mat'); 
end

[C.LMcells,C.V1cells,C.Behcells,C.params] = cleanupcells(cleanupcriteria,C.LMcells,C.V1cells,C.Behcells,C.params);
[cell_keep_ind,C.LMcells,C.V1cells,C.Behcells,C.params] = removecells(cell_keep_ind,cellremovecriteria,C.LMcells,C.V1cells,C.Behcells,C.params);

[IC.LMcells,IC.V1cells,IC.Behcells,IC.params] = cleanupcells(cleanupcriteria,IC.LMcells,IC.V1cells,IC.Behcells,IC.params);
[cell_keep_ind,IC.LMcells,IC.V1cells,IC.Behcells,params] = removecells(cell_keep_ind,cellremovecriteria,IC.LMcells,IC.V1cells,IC.Behcells,IC.params);



prctmax = 0.1;
visualize_spike_classes = 1;
[C.LMcells,C.V1cells,C.Behcells,C.params] = get_spikeshape_properties(prctmax,C.LMcells,C.V1cells,C.Behcells,C.params);
[IC.LMcells,IC.V1cells,IC.Behcells,IC.params] = get_spikeshape_properties(prctmax,IC.LMcells,IC.V1cells,IC.Behcells,IC.params);


if strcmp(exptype,'FF')
    targetcell = C.LMcells;
    ICtargetcells = IC.LMcells;
else
    targetcell = C.V1cells;
    ICtargetcells = IC.V1cells;
end
%% IMPORTANT:to run this, load cells (def 150ms),  assign spikeshapes, set targetcell below, then run
% clean and remove according to the criteria: default: R0C1
include_phys = 1;

% smb or smb_ls
allgo = cell2mat(cellfun(@(x) x.smb_ls.go, targetcell,'UniformOutput',0));
allnogo = cell2mat(cellfun(@(x) x.smb_ls.nogo, targetcell,'UniformOutput',0));

allgobs = cell2mat(cellfun(@(x) x.smb_bs.go, targetcell,'UniformOutput',0));
allnogobs = cell2mat(cellfun(@(x) x.smb_bs.nogo, targetcell,'UniformOutput',0));

if include_phys
    allretinx=cell2mat(cellfun(@(y) repmat(y,1,8),cellfun(@(x) {x.retin.center(1)}, targetcell),'UniformOutput',0));
    allretiny=cell2mat(cellfun(@(y) repmat(y,1,8),cellfun(@(x) {x.retin.center(1)}, targetcell),'UniformOutput',0));
    % 2 options: can change change from .distancelocal (distance to center of local shanks) to .distancelong
    % (distance to center of shanks in the other area)
    allretinlocald=cell2mat(cellfun(@(y) repmat(y,1,8),cellfun(@(x) {x.retin.distancelong}, targetcell),'UniformOutput',0));
    
    allAvfr=cell2mat(cellfun(@(y) repmat(y,1,8),cellfun(@(x) {x.Avfr},targetcell),'UniformOutput',0));
    alldepth=cell2mat(cellfun(@(y) repmat(y,1,8),cellfun(@(x) {x.depth},targetcell),'UniformOutput',0));
    
    % first change this, make sure others do or dont
    allspikeamp=cell2mat(cellfun(@(y) repmat(y,1,8),cellfun(@(x) {x.spikeamplitude},targetcell),'UniformOutput',0));
    allspikewidth=cell2mat(cellfun(@(y) repmat(y,1,8),cellfun(@(x) {x.spikewidth},targetcell),'UniformOutput',0));
    % index of CENTERS
    alltimecat = cell2mat(cellfun(@(x) {[1:8]},targetcell));
    
    
    monitormiddlex = 3;
    monitormiddley = 2.5;
    allcenterdistance = sqrt( (allretinx - monitormiddlex).^2 + (allretiny - monitormiddley).^2 );
end

%% behavioral params
% lick:
% figure;scatter(cellfun(@(x,y) SI(mean(x.stimspikes.go),mean(y.stimspikes.go)),targetcell,ICtargetcells),cellfun(@(x,y) SI(mean(y.stimspikes.nogo),mean(x.stimspikes.nogo)),targetcell,ICtargetcells))
notenough =[];%find(cellfun(@(x) numel(x.stimspikes.nogo)<5 ||  numel(x.stimspikes.go)<5,ICtargetcells));
% averaged for go and nogo
licktun_percell = nanmean([cellfun(@(x,y) SI(mean(x.stimspikes.go),mean(y.stimspikes.go)),targetcell,ICtargetcells);cellfun(@(x,y) SI(mean(y.stimspikes.nogo),mean(x.stimspikes.nogo)),targetcell,ICtargetcells)],1);
%licktun_percell = cellfun(@(x,y) SI(mean(x.stimspikes.go),mean(y.stimspikes.go)),targetcell,ICtargetcells);
licktun_percell(notenough) = nan;
all_lick_tun= cell2mat(cellfun(@(y) repmat(y,1,8),arrayfun(@(x) {x}, licktun_percell),'UniformOutput',0));
%figure;histogram(licktun_percell,[-1:0.1:1],'Normalization','pdf'); ylabel('lick tuning')

earlystimtun_percell = abs(cellfun(@(x) SI(nanmean(x.laAbs.go{1}),nanmean(x.laAbs.nogo{1})),targetcell));
%figure;hist(earlystimtun_percell)
earlystimtun_percell(notenough) = nan;
all_earlystim_tun=cell2mat(cellfun(@(y) repmat(y,1,8),arrayfun(@(x) {x}, earlystimtun_percell),'UniformOutput',0));


earlystimtun_percell_sign = sign(cellfun(@(x) SI(nanmean(x.laAbs.go{1}),nanmean(x.laAbs.nogo{1})),targetcell));
earlystimtun_percell_sign(notenough) = nan;
%figure;hist(earlystimtun_percell)
all_earlystim_tun_sign=cell2mat(cellfun(@(y) repmat(y,1,8),arrayfun(@(x) {x}, earlystimtun_percell_sign),'UniformOutput',0));


stimtun_percell = abs(nanmean([cellfun(@(x,y) SI(mean(x.stimspikes.go),mean(y.stimspikes.nogo)),targetcell,ICtargetcells);cellfun(@(x,y) SI(mean(y.stimspikes.go),mean(x.stimspikes.nogo)),targetcell,ICtargetcells)],1));
stimtun_percell(notenough) = nan;
all_stim_tun= cell2mat(cellfun(@(y) repmat(y,1,8),arrayfun(@(x) {x}, stimtun_percell),'UniformOutput',0));

stimtun_percell_sign = sign(nanmean([cellfun(@(x,y) SI(mean(x.stimspikes.go),mean(y.stimspikes.nogo)),targetcell,ICtargetcells);cellfun(@(x,y) SI(mean(y.stimspikes.go),mean(x.stimspikes.nogo)),targetcell,ICtargetcells)],1));
stimtun_percell_sign(notenough) = nan;
all_stim_tun_sign= cell2mat(cellfun(@(y) repmat(y,1,8),arrayfun(@(x) {x}, stimtun_percell_sign),'UniformOutput',0));



%%%%% savibg vars for histograms: re is pre removing nans universally
if strcmp(exptype,'FF')
    FF.pre.all_lick_tun = all_lick_tun(find(~isnan(all_lick_tun)));
    FF.pre.all_earlystim_tun = all_earlystim_tun(find(~isnan(all_earlystim_tun)));
    FF.pre.all_earlystim_tun_sign = all_earlystim_tun_sign(find(~isnan(all_earlystim_tun_sign)));
    FF.pre.all_stim_tun = all_stim_tun(find(~isnan(all_stim_tun)));
    FF.pre.all_stim_tun_sign = all_stim_tun_sign(find(~isnan(all_stim_tun_sign)));
    
else
    FB.pre.all_lick_tun = all_lick_tun(find(~isnan(all_lick_tun)));
    FB.pre.all_earlystim_tun = all_earlystim_tun(find(~isnan(all_earlystim_tun)));
    FB.pre.all_earlystim_tun_sign = all_earlystim_tun_sign(find(~isnan(all_earlystim_tun_sign)));
    FB.pre.all_stim_tun = all_stim_tun(find(~isnan(all_stim_tun)));
    FB.pre.all_stim_tun_sign = all_stim_tun_sign(find(~isnan(all_stim_tun_sign)));
end


%% remove all nans and infs
if include_phys
    nonanind = find(~isnan(mean([allretinlocald;allcenterdistance;allAvfr;alldepth;allspikeamp;allspikewidth;all_lick_tun;all_earlystim_tun;all_earlystim_tun_sign;all_stim_tun;all_stim_tun_sign],1)));
    allretinlocald = allretinlocald(nonanind);
    allcenterdistance = allcenterdistance(nonanind);
    allAvfr = allAvfr(nonanind);
    alldepth = alldepth(nonanind);
    allspikeamp = allspikeamp(nonanind);
    allspikewidth = allspikewidth(nonanind);
    alltimecat = alltimecat(nonanind);
    all_lick_tun = all_lick_tun(nonanind);
    all_earlystim_tun = all_earlystim_tun(nonanind);
    all_earlystim_tun_sign = all_earlystim_tun_sign(nonanind);
    all_stim_tun = all_stim_tun(nonanind);
    all_stim_tun_sign = all_stim_tun_sign(nonanind);
else
    nonanind = find(~isnan(mean([all_lick_tun;all_earlystim_tun;all_earlystim_tun_sign;all_stim_tun;all_stim_tun_sign],1)));
    all_lick_tun = all_lick_tun(nonanind);
    all_earlystim_tun = all_earlystim_tun(nonanind);
    all_earlystim_tun_sign = all_earlystim_tun_sign(nonanind);
    all_stim_tun = all_stim_tun(nonanind);
    all_stim_tun_sign = all_stim_tun_sign(nonanind);
end


allgo = allgo(nonanind);
allnogo = allnogo(nonanind);
allgobs = allgobs(nonanind);
allnogobs = allnogobs(nonanind);
%%% save arrays after removing nans
if strcmp(exptype,'FF')
    FF.post.all_lick_tun = all_lick_tun;
    FF.post.all_earlystim_tun = all_earlystim_tun;
    FF.post.all_earlystim_tun_sign = all_earlystim_tun_sign;
    FF.post.all_stim_tun = all_stim_tun;
    FF.post.all_stim_tun_sign = all_stim_tun_sign;

else
    FB.post.all_lick_tun = all_lick_tun;
    FB.post.all_earlystim_tun = all_earlystim_tun;
    FB.post.all_earlystim_tun_sign = all_earlystim_tun_sign;
    FB.post.all_stim_tun = all_stim_tun;
    FB.post.all_stim_tun_sign = all_stim_tun_sign;
end


%% Example exploratory plots, relation of receptive field location and silencing

% figure;scatter(allretinx,allretiny,100,allnogo,'.');xlim([1,5]);ylim([1,4])
% figure;scatter([allretinlocald,allretinlocald],[allgo,allnogo],100,'k.')
% figure;subplot(1,2,1);scatter([allretinx,allretinx],[allgo,allnogo],100,'k.')
% subplot(1,2,2);scatter([allgo,allnogo],[allretiny,allretiny],100,'k.')
%%

rngnum = 5; 
nfold = 20; 
norm2obv = 1; 
% Keep 1 in order to be able to compare.

glmprop.link = 'identity';
glmprop.spec = 'quadratic';
glmprop.cat = []; % index of categorical predictors: default: []

figure;

plotprop.plotprediction = 0;
plotprop.plotLL = 1;
plotprop.LLhandle = gca;
plotprop.LLloc = 1;
plotprop.singlectrl = 0; 
plotprop.scatterpl = 0;
plotprop.comparisons = 'perfold'; % 'perfold' or 'perobsv': discribed in methods. Perobs: not sure if complemety correct

tasktype = 'all'; % all, go, nogo


%%% prepare matrices:
clear X Y
Y.go = allgo; Y.nogo = allnogo;

X.go = allgobs; X.nogo = allnogobs;
plotprop.LLloc = 1;
[outmodels,Normalized_logl1] = predict_Slc_static(X,Y,tasktype,rngnum,nfold,norm2obv,glmprop,plotprop);
if include_phys
    %%%%%%%%%%%%%%%% retin local distance -- or long range distance (depending on the first section)
    X.go = [allgobs;allretinlocald]; X.nogo = [allnogobs;allretinlocald];
    plotprop.LLloc = 2;
    [outmodels,Normalized_logl2] = predict_Slc_static(X,Y,tasktype,rngnum,nfold,norm2obv,glmprop,plotprop);
    
    if plotprop.singlectrl
        X.go = allretinlocald; X.nogo = allretinlocald;
        plotprop.LLloc =2.2;
        [outmodels,Normalized_logl] = predict_Slc_static(X,Y,tasktype,rngnum,nfold,norm2obv,glmprop,plotprop);
    end
    %%%%%%%%%%%%%%% retin center monitor distance
    X.go = [allgobs;allcenterdistance]; X.nogo = [allnogobs;allcenterdistance];
    plotprop.LLloc = 3;
    [outmodels_r,Normalized_logl3] = predict_Slc_static(X,Y,tasktype,rngnum,nfold,norm2obv,glmprop,plotprop);
    
    if plotprop.singlectrl
        X.go = allcenterdistance; X.nogo = allcenterdistance;
        plotprop.LLloc =3.2;
        [outmodels,Normalized_logl] = predict_Slc_static(X,Y,tasktype,rngnum,nfold,norm2obv,glmprop,plotprop);
        
    end
    %%%%%%%%%%%%%%% Avfr
    X.go = [allgobs;allAvfr]; X.nogo = [allnogobs;allAvfr];
    plotprop.LLloc = 4;
    [outmodels,Normalized_logl4] = predict_Slc_static(X,Y,tasktype,rngnum,nfold,norm2obv,glmprop,plotprop);
    
    if plotprop.singlectrl
        X.go = allAvfr; X.nogo = allAvfr;
        plotprop.LLloc = 4.2;
        [outmodels,Normalized_logl] = predict_Slc_static(X,Y,tasktype,rngnum,nfold,norm2obv,glmprop,plotprop);
    end
    %%%%%%%%%%%%%%% depth
    plotprop.plotprediction = 0;
    
    X.go = [allgobs;alldepth]; X.nogo = [allnogobs;alldepth];
    plotprop.LLloc = 5;
    [outmodels_d,Normalized_logl5] = predict_Slc_static(X,Y,tasktype,rngnum,nfold,norm2obv,glmprop,plotprop);
    
    if plotprop.singlectrl
        X.go = alldepth; X.nogo = alldepth;
        plotprop.LLloc =5.2;
        [outmodels,Normalized_logl] = predict_Slc_static(X,Y,tasktype,rngnum,nfold,norm2obv,glmprop,plotprop);
        
    end
    plotprop.plotprediction = 0;
    
    
    %%%%%%%%%%%%%%% spike width
    X.go = [allgobs;allspikewidth]; X.nogo = [allnogobs;allspikewidth];
    plotprop.LLloc = 6;
    [outmodels,Normalized_logl6] = predict_Slc_static(X,Y,tasktype,rngnum,nfold,norm2obv,glmprop,plotprop);
    
    if plotprop.singlectrl
        X.go = allspikewidth; X.nogo = allspikewidth;
        plotprop.LLloc =6.2;
        [outmodels,Normalized_logl] = predict_Slc_static(X,Y,tasktype,rngnum,nfold,norm2obv,glmprop,plotprop);
        
    end
    %%%%%%%%%%%%%%% spike amplitude
    X.go = [allgobs;allspikeamp]; X.nogo = [allnogobs;allspikeamp];
    plotprop.LLloc = 7;
    [outmodels,Normalized_logl] = predict_Slc_static(X,Y,tasktype,rngnum,nfold,norm2obv,glmprop,plotprop);
    
    if plotprop.singlectrl
        X.go = allspikeamp; X.nogo = allspikeamp;
        plotprop.LLloc =7.2;
        [outmodels,Normalized_logl] = predict_Slc_static(X,Y,tasktype,rngnum,nfold,norm2obv,glmprop,plotprop);
        
    end
    
    
    %%%%%%%%%%%%%%% categoical: tasktype: only when all
    glmprop.cat = [];
    X.go = [allgobs;zeros(size(allgobs))];X.nogo = [allnogobs;ones(size(allnogobs))];
    plotprop.LLloc = 8;
    [outmodels,Normalized_logl8] = predict_Slc_static(X,Y,tasktype,rngnum,nfold,norm2obv,glmprop,plotprop);
    glmprop.cat = [];
    
    if plotprop.singlectrl
        glmprop.cat = 1;
        X.go = zeros(size(allgobs)); X.nogo = zeros(size(allnogobs));
        plotprop.LLloc = 8.2;
        [poutmodels,Normalized_logl] = predict_Slc_static(X,Y,tasktype,rngnum,nfold,norm2obv,glmprop,plotprop);
        glmprop.cat = [];
    end
    
    
    %%%%%%%%%%%%%% catergorical: time
    plotprop.plotprediction = 0;
    glmprop.cat = [];
    X.go = [allgobs;alltimecat]; X.nogo = [allnogobs;alltimecat];
    plotprop.LLloc = 9;
    [outmodels,Normalized_logl9] = predict_Slc_static(X,Y,tasktype,rngnum,nfold,norm2obv,glmprop,plotprop);
    %glmprop.link = 'identity';
    glmprop.cat = [];
    
    if plotprop.singlectrl
        glmprop.cat = [];
        X.go = [alltimecat]; X.nogo = [alltimecat];
        plotprop.LLloc = 9.2;
        [outmodels,Normalized_logl] = predict_Slc_static(X,Y,tasktype,rngnum,nfold,norm2obv,glmprop,plotprop);
        glmprop.cat = [];
    end
    plotprop.plotprediction = 0;
end
%%%%%%%%%%%%%%% licktuning
X.go = [allgobs;all_lick_tun]; X.nogo = [allnogobs;all_lick_tun];
plotprop.LLloc = 10;
[outmodels,Normalized_logl10] = predict_Slc_static(X,Y,tasktype,rngnum,nfold,norm2obv,glmprop,plotprop);

if plotprop.singlectrl
    X.go = all_lick_tun; X.nogo = all_lick_tun;
    plotprop.LLloc =10.2;
    [outmodels,Normalized_logl] = predict_Slc_static(X,Y,tasktype,rngnum,nfold,norm2obv,glmprop,plotprop);
    
end
%%%%%%%%%%%%%%% earlystimtuning
X.go = [allgobs;all_earlystim_tun]; X.nogo = [allnogobs;all_earlystim_tun];
plotprop.LLloc = 11;
[outmodels,Normalized_logl11] = predict_Slc_static(X,Y,tasktype,rngnum,nfold,norm2obv,glmprop,plotprop);

if plotprop.singlectrl
    X.go = all_earlystim_tun; X.nogo = all_earlystim_tun;
    plotprop.LLloc =11.2;
    [outmodels,Normalized_logl] = predict_Slc_static(X,Y,tasktype,rngnum,nfold,norm2obv,glmprop,plotprop);
    
end
%%%%%%%%%%%%%%% earlystimtuning: categorical

glmprop.cat = [2];
X.go = [allgobs;all_earlystim_tun_sign]; X.nogo = [allnogobs;all_earlystim_tun_sign];
plotprop.LLloc = 12;
[outmodels,Normalized_logl9] = predict_Slc_static(X,Y,tasktype,rngnum,nfold,norm2obv,glmprop,plotprop);
%glmprop.link = 'identity';
glmprop.cat = [];

if plotprop.singlectrl
    glmprop.cat = [];
    X.go = [all_earlystim_tun_sign]; X.nogo = [all_earlystim_tun_sign];
    plotprop.LLloc = 12.2;
    [outmodels,Normalized_logl] = predict_Slc_static(X,Y,tasktype,rngnum,nfold,norm2obv,glmprop,plotprop);
    glmprop.cat = [];
end

%%%%%%%%%%%%%%% stimtuning
X.go = [allgobs;all_stim_tun]; X.nogo = [allnogobs;all_stim_tun];
plotprop.LLloc = 13;
[outmodels,Normalized_logl13] = predict_Slc_static(X,Y,tasktype,rngnum,nfold,norm2obv,glmprop,plotprop);

if plotprop.singlectrl
    X.go = all_stim_tun; X.nogo = all_stim_tun;
    plotprop.LLloc =13.2;
    [outmodels,Normalized_logl] = predict_Slc_static(X,Y,tasktype,rngnum,nfold,norm2obv,glmprop,plotprop);
    
end

%%%%%%%%%%%%%%% stimtuning: categorical
glmprop.cat = [2];
X.go = [allgobs;all_stim_tun_sign]; X.nogo = [allnogobs;all_stim_tun_sign];
plotprop.LLloc = 14;
[outmodels,Normalized_logl14] = predict_Slc_static(X,Y,tasktype,rngnum,nfold,norm2obv,glmprop,plotprop);
glmprop.cat = [];

if plotprop.singlectrl
    X.go = all_stim_tu_signn; X.nogo = all_stim_tun_sign;
    plotprop.LLloc =13.2;
    [outmodels,Normalized_logl] = predict_Slc_static(X,Y,tasktype,rngnum,nfold,norm2obv,glmprop,plotprop);
    
end

    

%%%%%%%%%%%
ylim(plotprop.LLhandle,[-0.5 0.5])% 0.05
xlim(plotprop.LLhandle,[0 plotprop.LLloc+1])


set(plotprop.LLhandle,'XTick',[1 2 2.1 3 3.1 4 4.1 5 5.1 6 6.1 7 7.1 8 8.1 9 9.1  10 10.1  11 11.1 12 12.1 13 13.1 14 14.1])
set(plotprop.LLhandle,'XTickLabel',{'bs','longDist+bs','chance','centerDist+bs','chance',...
    'Avfr+bs','chance','depth+bs','chance',...
    'spikeWidth+bs','chance', 'spikeAmp+bs','chance',...
    'task+bs','chance','time+bs','chance',...
    'lick-tun+bs','chance','earlystim-tun+bs','chance','earlystimid-tun+bs','chance','stim-tun+bs','chance','stim-id+bs','chance'})
set(plotprop.LLhandle,'XTickLabelRotation',90)

%%%%
