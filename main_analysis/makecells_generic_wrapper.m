clear all
%%
addpath(genpath('/mnt/data/Mitra/cache/codes/matlabkwikreader/'))
addpath(genpath('/mnt/data/Mitra/cache/codes/'))
%%
% Animals are selected based on performance. d-prime>threshold. Mice with d
% prime lower than the threshold are excluded from any further analysis.

% animallist = {'MPV17','MPV18_2',...
%     'VL53','VL52','VL51','VL66','MPV35_2'};%,...
%      %'MPV36_2','MPV30'};    
% preprocessinglist = {'2020_03_02_13_10_37','2020_03_02_13_27_45',...
%     '2017_12_05_18_58_17','2017_12_04_15_08_47','2017_12_03_17_30_15','2018_03_20_11_53_19','2020_03_02_14_31_53',...
%     };
% 
% exptype = {'FB','FB',...
%     'FB','FB','FB','FB','FB'};
% prefix = '';


animallist ={'VL61','VL63','VL55','VL59',...
    'MPV33','MPV31','MPV34_2'};
preprocessinglist = {'2018_03_16_16_53_13','2018_04_01_12_24_38','2018_04_01_15_03_48','2018_04_01_13_24_48',...
    '2020_03_02_19_35_12','2020_03_02_19_58_50','2020_03_02_20_32_35'};

exptype = {'FF','FF','FF','FF',...
    'FF','FF','FF'};
prefix = '';

%% makecells parameters
% IMPORTANT: Window is assumed symmetric for most purposes : half of PAllOn is
% considered point zero

clear params

params.samplingf = 30000;
params.Window = 30000;


params.min_n_trials_per_delay = 15;%10 when excludegrooming=0. set to 15 for FF
%
params.excludegrooming = 1; %def: 1 -1:only grooming. When this is -1, only correct should be 0. For correct trials, it is 1
params.onlycorrect = 1; % def: 1 change to 1 for correct trials only.
params.onlystableperiod = 0; %def:0

% (class 2 params) plotting
params.edgestep_pl = 20*30;%
params.alignstyle_pl = 1;

% (class 3 params) analysis
params.edgestep_an = 80*30;%80*30 or 150*30
params.onsetlag_an = nan;%40*30;
% (class 4 params) parameters for moving window with lag from laser onset (effect onset
% detection)
params.movingwin_withlags.lagsize = 5*30;% 5
params.movingwin_withlags.nlags = 0; %40
params.movingwin_withlags.movingwinsize = 20*30; % 20*30
%% saving parameters
saveoutput = 1;
savedir = '/mnt/data/Mitra/figs/P2_L/bothdircombined';
% change name if FF and FB combined
%savename = sprintf('cells_FB_%dms_lw%dms_exG%d_onlyC%d_onlyS%d.mat',params.edgestep/30,params.movingwin_withlags.movingwinsize/30,params.excludegrooming,params.onlycorrect,params.onlystableperiod);
savename = sprintf('cells%s_%s_pl%d_an%d_lw%d_exG%d_onlyC%d_onlyS%d_plstyle%d.mat',...
    prefix,exptype{1},params.edgestep_pl/30,params.edgestep_an/30,params.movingwin_withlags.movingwinsize/30,...
    params.excludegrooming,params.onlycorrect,params.onlystableperiod,...
    params.alignstyle_pl);
if isfield(params,'onsetlag_an') && ~isnan(params.onsetlag_an)
    savename = sprintf('cells%s_%s_pl%d_an%d-%d_lw%d_exG%d_onlyC%d_onlyS%d_plstyle%d.mat',...
    prefix,exptype{1},params.edgestep_pl/30,params.onsetlag_an/30,params.edgestep_an/30,...
    params.movingwin_withlags.movingwinsize/30,...
    params.excludegrooming,params.onlycorrect,params.onlystableperiod,...
    params.alignstyle_pl);
end

%%
trialtypelist = {'gotrialind','nogotrialind'};
% selecting target trials based on params
if params.onlycorrect == 1
    targettrialtypelist = {'correctgotrialind','correctnogotrialind'};
elseif params.onlycorrect == -1
    targettrialtypelist = {'incorrectgotrialind','incorrectnogotrialind'};
else
    targettrialtypelist = trialtypelist;
end
tic
[LMcells,V1cells,Behcells] = makecells_generic_core(params,animallist,preprocessinglist,exptype,trialtypelist,targettrialtypelist);
toc
if saveoutput
    cd(savedir)
    save(savename,'LMcells','V1cells','params','Behcells')
end
