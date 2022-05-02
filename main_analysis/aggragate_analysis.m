% load cells_xms
%%
addpath(genpath('/mnt/data/Mitra/cache/codes/'))
rmpath(genpath('/mnt/data/Mitra/cache/codes/chronux_2_12/')); % has annoying funcs, similar names as matlab natives
cd('/mnt/data/Mitra/figs/P2_L/bothdircombined')
%% global plot parameters
globalplot.histnorm = 'pdf'; 

%% choose params for loading cells
% load feedback or feedforward
ftl.savenameprefix = 'FF'; % 'FB' or 'FF'
ftl.taskprefix = ''; %or {'_U'}
ftl.pl_ms = 40; % plotting bin size
ftl.an_ms = 150;
ftl.lw_ms = 20;
ftl.exG = 1;
ftl.onlyC = 1;
ftl.onlyS = 0;
ftl.plstyle = 1; % plot style: 0 is old
cellfiletoload = sprintf('cells%s_%s_pl%d_an%d_lw%d_exG%d_onlyC%d_onlyS%d_plstyle%d.mat',...
    ftl.taskprefix,ftl.savenameprefix,...
    ftl.pl_ms,ftl.an_ms,ftl.lw_ms,...
    ftl.exG,ftl.onlyC,ftl.onlyS,...
    ftl.plstyle);

%% load, clean up and cell remove params

cleanupcriteria = struct;
cleanupcriteria.doclean = 1; % 1 or 0, perform any cleaning or not. If zero does nothing
cleanupcriteria.significance = nan; % significance of activity instead of absolute value
cleanupcriteria.smb_activity = 2.5; % in Hz

%%%%
%%% TODO: play with responsivenes methods and criteria, look at plots and
%%% improve methods and params choice. 3 is tricky: removed late ramping
%%% high baseline neurons

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

%%% clean and remove cells


[LMcells,V1cells,Behcells,params] = cleanupcells(cleanupcriteria,LMcells,V1cells,Behcells,params);
[cell_keep_ind,LMcells,V1cells,Behcells,params] = removecells(cell_keep_ind,cellremovecriteria,LMcells,V1cells,Behcells,params);


%% cells are finalized
%%%%% from this point on, the number of cells is the same - subselecting
%%%%% doesn't remove cells
%% check equal number of laser centers: validity of the loaded cell structure
[V1cells,LMcells] = SubselectExpType(V1cells,LMcells,ftl.savenameprefix);

%% check spike shapes: get width: static property of cells: how does it interact with cell indexing? add to fields of the cell?
% assign width and amplitude to all cells
prctmax = 0.1;
visualize_spike_classes = 1;
[LMcells,V1cells,Behcells,params] = get_spikeshape_properties(prctmax,LMcells,V1cells,Behcells,params);

% to test prctmax, and set th and th2, you can play with params and
% visualize spikes in this script
% just visualization, V1cells, LMcells dont change
th=12;%9
th2 = 18;%18
if visualize_spike_classes
    spikeshape_classification
end

%% static measures For both FF and FB: agnostic to FF and FB, just a choice of V1 or LM cells
% static analysis 1 - onset and offset of silencing
%%%%%%%%%%%%%%% set params
% lw bin is the only important bit for calculations. 
% cellfiletoload = 'cells_FB_40ms_lw20ms_exG1_onlyC1_onlyS0';
% choose V1 or LM, and the effect you want to look at
targetcellname = 'LMcells';
effect = 'onset'; % 'onset' or 'offset'
subselectcells = 0;
showspikeshape = 1;
% selecting based on depth
%subselect_ind_phrase = sprintf('find(cellfun(@(x) x.depth,%s)<400 & cellfun(@(x) x.depth,%s)>200)',targetcellname,targetcellname);
% spikeshape type 1
%subselect_ind_phrase = sprintf('find(cellfun(@(x) x.spikewidth,%s)>th & cellfun(@(x) x.spikewidth,%s)<th2)',targetcellname,targetcellname);
% spike shape type 2
%subselect_ind_phrase = sprintf('find(cellfun(@(x) x.spikewidth,%s)<=th)',targetcellname);

%%%%%%%%%%%%%%%

get_effect_onset;
% hists are for each cell, averaged over all silencing lags (calculation
% method is biased, can't be compared for different lags)
%%%%%%%%%%%%%%%%



% static analysis 2 - histogram of effects in both areas
% cellfiletoload = 'cells_FB_150ms_exG1_onlyC1_onlyS0.mat';
targetcellname = 'V1cells';
subselectcells = 0;
showspikeshape = 1;
% selecting based on depth
%subselect_ind_phrase = sprintf('find(cellfun(@(x) x.depth,%s)<400 & cellfun(@(x) x.depth,%s)>200)',targetcellname,targetcellname);
% spikeshape type 1
subselect_ind_phrase = sprintf('find(cellfun(@(x) x.spikewidth,%s)>th & cellfun(@(x) x.spikewidth,%s)<th2)',targetcellname,targetcellname);
% spike shape type 2
%subselect_ind_phrase = sprintf('find(cellfun(@(x) x.spikewidth,%s)<=th)',targetcellname);

static_effectsize_histograms;

%% For combining FF and FB run this, then select responsive cells, then plot traces
FFIC =load('cells_FF_pl10_an80_lw20_exG0_onlyC-1_onlyS0_plstyle1.mat');
FBIC =load('cells_FB_pl10_an80_lw20_exG0_onlyC-1_onlyS0_plstyle1.mat');
FFC =load('cells_FF_pl10_an80_lw20_exG1_onlyC1_onlyS0_plstyle1.mat');
FBC =load('cells_FB_pl10_an80_lw20_exG1_onlyC1_onlyS0_plstyle1.mat');
% %
V1cells = [FFC.V1cells, FBC.V1cells];
LMcells = [FFC.LMcells, FBC.LMcells];
params = FBC.params;
Behcells = FBC.Behcells;
% get cell_keep_ind according to correct trials
cell_keep_ind =nan;
[cell_keep_ind,LMcells,V1cells,Behcells,params] = ...
    removecells(cell_keep_ind,cellremovecriteria,LMcells,V1cells,Behcells,params);

% apply cell_keep_ind to incorrect trials
V1cells = [FFIC.V1cells, FBIC.V1cells];
LMcells = [FFIC.LMcells, FBIC.LMcells]; 
params = FBIC.params;
Behcells = FBIC.Behcells;
[cell_keep_ind,LMcells,V1cells,Behcells,params] = ...
    removecells(cell_keep_ind,cellremovecriteria,LMcells,V1cells,Behcells,params);

%% basline imagesc and trace plots: activity, zscore selectivity
plotstodo.trace = 1;
plotstodo.raster = 0;
plotstodo.sel = 0;
ntrial_th = 1; % keep 1., increase if you want to exclude animals with less trials
Aggregate_BaselineGovsNogoPlots(V1cells,LMcells,params,plotstodo,ntrial_th)

