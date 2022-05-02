clear all

%%
% Large vectors like PdTrig and LaserTrig are not saved. They can easily be
% recalculated later if needed.
%
%
%% default plot properties
set(groot,'defaultAxesTitleFontWeight','normal')
set(groot,'defaultAxesFontSize',6)

%%
addpath(genpath('/mnt/data/Mitra/cache/codes/matlabkwikreader/'))
addpath(genpath('/mnt/data/Mitra/cache/codes/'))
rmpath('/mnt/data/Mitra/cache/codes/chronux_2_12/spectral_analysis/continuous'); % confusing findpeaks in chronux
                                                                                 % that interferes with matlab's native one

%% params
make_params_list_MPV17
%% make saving folder
cd(savingdirroot)
% make folders if necessary
if ~exist(fullfile(savingdirroot,'preprocessing'),'dir')
    mkdir('preprocessing')
end
% timestamped name
cd('preprocessing')
savename=datestr(now,'yyyy_mm_dd_HH_MM_SS');
mkdir(savename);
%% constructing paths and assigning channels: be careful that it is /processed or /manual_sorting
% channels
Pdch = n_ephys_channles + n_aux_channels + photodiode_channel;
Laserch = n_ephys_channles + n_aux_channels + laser_channel;

% paths
kwdrecording = ['/recordings/',rec_number,'/data'];
%raw files:

path2kwe = ['/mnt/data/Mitra/raw/',animalname,expname{1},'/experiment1.kwe']; % for messages
path2kwd = ['/mnt/data/Mitra/raw/',animalname,expname{1},'/experiment1_',processor_number,'.raw.kwd']; %for pd % 101 usually

% spikes : first experiment in the list only
V1{1}.path2spikes = ['/mnt/data/Mitra/processed/',animalname,expname{1},'_kwd_',processor_number,'_rec',rec_number,'_shkV1_shk0/']; % kilosorted spikes
V1{2}.path2spikes =  ['/mnt/data/Mitra/processed/',animalname,expname{1},'_kwd_',processor_number,'_rec',rec_number,'_shkV1_shk1/'];
LM{1}.path2spikes =  ['/mnt/data/Mitra/processed/',animalname,expname{1},'_kwd_',processor_number,'_rec',rec_number,'_shkLM_shk0/'];
LM{2}.path2spikes =  ['/mnt/data/Mitra/processed/',animalname,expname{1},'_kwd_',processor_number,'_rec',rec_number,'_shkLM_shk1/'];
%% get info
[kwe_infos, messages, ttl] = read_kwe_file(path2kwe);
disp(messages(find(messages.nodeID==100),:).text);
start_time=messages(find(messages.nodeID==100),:).time_samples
%messages(find(messages.eventID==0),:).text
info = h5info(path2kwd,kwdrecording);
rec_length=info.Dataspace.Size(2);


%% get all spike times
if strcmp(units_to_extract,'gmu')
    [V1,LM] = Extractspikeinfo2_gmuonly(V1,LM);
elseif strcmp(units_to_extract,'all')
    [V1,LM] = Extractspikeinfo2(V1,LM);
end

V1_Units = struct;
LM_Units = struct;

%%% separate a subunit from MU based on depth here:
LM_Units.MUSpikes_AllExps=cell2mat([LM{1,1}.AllUnitSpikeTimes,LM{1,2}.AllUnitSpikeTimes]');
V1_Units.MUSpikes_AllExps=cell2mat([V1{1,1}.AllUnitSpikeTimes,V1{1,2}.AllUnitSpikeTimes]');

LM_Units.SingleUnitSpikeTimes_AllExps = [LM{1,1}.SingleUnitSpikeTimes, LM{1,2}.SingleUnitSpikeTimes];
V1_Units.SingleUnitSpikeTimes_AllExps = [V1{1,1}.SingleUnitSpikeTimes, V1{1,2}.SingleUnitSpikeTimes];

%% cutting all spikes to different experiments
end_points = rec_length;
start_points = 0 ;
%%% MU

LM_Units.MUSpikes = LM_Units.MUSpikes_AllExps(find((LM_Units.MUSpikes_AllExps<end_points) )) - start_points;
V1_Units.MUSpikes = V1_Units.MUSpikes_AllExps(find((V1_Units.MUSpikes_AllExps<end_points) )) - start_points;


%%% SU
for i = 1:length(V1_Units.SingleUnitSpikeTimes_AllExps)
    V1_Units.SUSpikes{i}= V1_Units.SingleUnitSpikeTimes_AllExps{i}(find((V1_Units.SingleUnitSpikeTimes_AllExps{i}<end_points) )) - start_points;
end

for i = 1:length(LM_Units.SingleUnitSpikeTimes_AllExps)
    LM_Units.SUSpikes{i}= LM_Units.SingleUnitSpikeTimes_AllExps{i}(find((LM_Units.SingleUnitSpikeTimes_AllExps{i}<end_points) )) - start_points;
end
%%% add spike shape, autocorr,...

%% parse ssmessages

testofflinemessages_MPV

sprintf('%d go trials from trial %d to %d',size(gotrialind,1),gotrialind(1),gotrialind(end))
sprintf('%d nogo trials from trial %d to %d',size(nogotrialind,1),nogotrialind(1),nogotrialind(end))

stabletrialind = 1:max(gotrialind(end),nogotrialind(end));

%% get pd andlaser trigger times
%%%%%%%%%%%%%%%%%% Determine start and end times, if not done already:
if isnan(tStart) |  isnan(tEnd)
    PdTrig_nofilter = h5read(path2kwd,kwdrecording,[Pdch 1],[1 inf]);
    figure;plot(PdTrig_nofilter(1:2*10^6));
    userin = input('tStart:','s');
    tStart = str2num(userin)
    
    figure;plot(PdTrig_nofilter(end-10^6:end));
    userin = input('How much to cut:','s');
    tEnd = length(PdTrig_nofilter) - str2num(userin)
end

%%%%%%%%%%%%%%%%%%%% get pd and laser onset/offset times

LaserTrig = h5read(path2kwd,kwdrecording,[Laserch 1],[1 inf]);
if strcmp(lasershape ,'ramp')
    [LStepTimeOn,LStepTimeStampOn] = FindEdges_rampdown(LaserTrig,DiffThL,tStart,tEnd);
elseif strcmp(lasershape,'square')
    [LStepTimeOn,~,LStepTimeStampOn,~] = FindEdges(LaserTrig,DiffThL,tStart,tEnd);
else
    disp('choose a lasershape first')
end
sprintf('%d laser onset detected',length(LStepTimeStampOn))
%  sprintf('%d laser offset detected',length(LStepTimeStampOff))

Window=WindowSec*samplingf;
LAllOn= repmat(LStepTimeStampOn',[1,2*Window+1])+repmat(double(-Window:Window),[length(LStepTimeStampOn'),1]); % if you don't double, it will be single and the subtraction would be fucked
pdfig_handle = figure('Units','normalized','Position',[0.2 0.1 0.7 0.7]);
subplot(3,1,1);plot((-Window:Window)/samplingf,LaserTrig(LAllOn(:,:))');
%%%
PdTrig_nofilter = h5read(path2kwd,kwdrecording,[Pdch 1],[1 inf]);
PdTrig=low_pass_filter(PdTrig_nofilter,samplingf,pdfiltfreq);
[PStepTimeOn,PStepTimeOff,PStepTimeStampOn,PStepTimeStampOff] = FindEdges(PdTrig,DiffTh,tStart,tEnd);

sprintf('%d pd onset detected',length(PStepTimeStampOn))
sprintf('%d pd offset detected',length(PStepTimeStampOff))

Window=WindowSec*samplingf;
WindowL=WindowSec*samplingf;
PAllOn= repmat(PStepTimeStampOn',[1,(Window+WindowL)+1])+repmat(double(-WindowL:Window),[length(PStepTimeStampOn'),1]); % if you don't double, it will be single and the subtraction would be fucked
subplot(3,1,2);plot((-WindowL:Window)/samplingf,PdTrig(PAllOn(:,:))');
subplot(3,1,3);plot((-WindowL:Window)/samplingf,LaserTrig(PAllOn(:,:))');
pdfig_handle.Children(3).Title.String = sprintf('%d pd onset detected\n%d pd offset detected\n%d laser onset detected',...
    length(PStepTimeStampOn),length(PStepTimeStampOff),length(LStepTimeStampOn));
pdfig_handle.Children(3).Title.FontWeight = 'normal';
pdfig_handle.Children(3).Title.FontSize = 9;

figure;plot((-WindowL:Window)/samplingf,PdTrig_nofilter(PAllOn(:,:))')


%% binning laser delays
[laserdelayfig_handle,LaserDelay,centers,LaserDelayBinned] = Extractlaserdelaynan(LStepTimeOn,PAllOn,Window,samplingf,Laserdelaybinwindow,Laserdelaybinnum);

%% behavior analysis
onlystableperiod = 0;
behavior_MPV;


%% check for stability (beh and units)
%%%% examples for intervals: [ 1:20 , 30:40]
%%%% 300:end   1:end  

fstab = figure;fstab.Units = 'normalized';fstab.Position = [0.1 0.1 0.8 0.8];
figure(fstab);ax=subplot(3,3,1);plot(movmean(correcttrials, [0 20]));ax.Title.String = 'perf over trials';
% MU
onlystableunits = 0;
checkmultiunitstability
% cut indices if necessary
prompt = 'want to cut?[y/n] ';
userin = input(prompt,'s');
if strcmp(userin,'y') || strcmp (userin,'Y')
    userin = input('enter stabletrialind in the form of indices:','s');
    allinds =  1:max(gotrialind(end),nogotrialind(end));
    sprintf('allinds(%s)',userin)
    stabletrialind = eval(sprintf('allinds(%s)',userin));
end

%% quality measure for each units for both cut session and whole session.
%%% adds a field to each area{shank} describing the unit drift over time. Both, through
%%% the whole experiment or in stabletrialind
onlystableperiod = 0;
checksingleunitstability

onlystableperiod = 1;
checksingleunitstability

% taking only high quality units, shows how  MU look over time
onlystableunits = 1;
driftthreshold = 0.2;
checkmultiunitstability;
%% in case of cutting, repeat behavior lick laser figures
%%% This is the last part of behavior_MPV, runs with the same behavior
%%% params otherwise, as before
onlystableperiod = 1;

normalizeforcorrect = 1;
figperfn_stab = figure; figperfn_stab.Name = 'normalized';
figperfn_stab.Units = 'normalized';
figperfn_stab.Position = [0.2 0.2 0.8 0.8];
makelaserlagperfplots;

normalizeforcorrect = 0;
figperf_stab = figure; figperf_stab.Name = 'not normalized';
figperf_stab.Units = 'normalized';
figperf_stab.Position = [0.2 0.2 0.8 0.8];
makelaserlagperfplots;

%% saving (pool from behavior)
prompt = 'Are you happy with preprocessing results? Y/N [Y]: ';
userin = input(prompt,'s');
if strcmp(userin,'y') || strcmp (userin,'Y')
    cd(savingdirroot)
    cd('preprocessing')
    cd(savename);
    % save variables
    save([expname{1}(2:end),'_preprocessing_',savename,'_ks2_full'],...
        'animalname','oemessagesfilename','V1','V1_Units','LM','LM_Units','units_to_extract',...
        'Window','WindowL','WindowSec','ttl','tStart','tEnd','start_time','start_points',...
        'savingdirroot','samplingf','rec_number','rec_length','path2kwd','path2kwe','expname','DiffTh','DiffThL',...
        'processor_number','photodiode_channel','Netevent_nodeID','n_ephys_channles','n_aux_channels','Pdch',...
        'Laserch','laser_channel','kwdrecording','kwe_infos','end_points',...
        'messages','LAllOn','PAllOn','pdfiltfreq','Laserdelaybinnum','Laserdelaybinwindow',...
        'LStepTimeOn','LStepTimeStampOn','centers','LaserDelay','LaserDelayBinned',...
        'PStepTimeOff','PStepTimeOn','PStepTimeStampOff','PStepTimeStampOn',...
        'SSmessages',...
        'gotrial','nogotrial','alltrialgn','gotrialind','nogotrialind','gostimid','nogostimid',...
        'lick_channel', 'running_channel','grooming_lick_threshold','nogroomingind',...
        'FA','Miss','correctnogotrialind','incorrectnogotrialind','correctgotrialind','incorrectgotrialind',...
        'speedWindowms', 'speedticks', 'speed',...
        'firstlicksample','firstvalvesample','lickcontforms','thr','valves','licks','lickthreshold','valvethreshold',...
        'allvalve','alllick','checkpos','lickthprctile','respwindowms','min_n_trials_per_delay',...
        'hitrate','missrate','crrate','farate','correctrate','incorrectrate',...
        'stabletrialind',...
        '-v7.3');
    
    %save figures
   % hgsave(pdfig_handle,[expname{1}(2:end),'_pd_',savename,'.fig'], '-v7.3');
    saveas(pdfig_handle,[expname{1}(2:end),'_pd_',savename,'.png']);
    
    hgsave(laserdelayfig_handle,[expname{1}(2:end),'_laserdelay_',savename,'.fig'], '-v7.3');
    saveas(laserdelayfig_handle,[expname{1}(2:end),'_laserdelay_',savename,'.png']);
    
    
    
    hgsave(beh_f1,[expname{1}(2:end),'_behavior_',savename,'.fig'], '-v7.3');
    saveas(beh_f1,[expname{1}(2:end),'_behavior_',savename,'.png']);
    
    hgsave(beh_f2,[expname{1}(2:end),'_behaviorexgrooming_',savename,'.fig'], '-v7.3');
    saveas(beh_f2,[expname{1}(2:end),'_behaviorexgrooming_',savename,'.png']);
    
    hgsave(f_lick1,[expname{1}(2:end),'_lick1_',savename,'.fig'], '-v7.3');
    saveas(f_lick1,[expname{1}(2:end),'_lick1_',savename,'.png']);
    
    hgsave(f_lick2,[expname{1}(2:end),'_lick2_',savename,'.fig'], '-v7.3');
    saveas(f_lick2,[expname{1}(2:end),'_lick2_',savename,'.png']);
    
    hgsave(f_behovertime,[expname{1}(2:end),'_perfmtime_',savename,'.fig'], '-v7.3');
    saveas(f_behovertime,[expname{1}(2:end),'_perfmtime_',savename,'.png']);
    
    hgsave(figperf,[expname{1}(2:end),'_laserperfmplots_',savename,'.fig'], '-v7.3');
    saveas(figperf,[expname{1}(2:end),'_laserperfmplots_',savename,'.png']);
    
    hgsave(figperfn,[expname{1}(2:end),'_laserperfmplotsnorm_',savename,'.fig'], '-v7.3');
    saveas(figperfn,[expname{1}(2:end),'_laserperfmplotsnorm_',savename,'.png']);
    
    hgsave(fstab,[expname{1}(2:end),'_spikebehaviorstab_',savename,'.fig'], '-v7.3');
    saveas(fstab,[expname{1}(2:end),'_spikebehaviorstab_',savename,'.png']);
    
    hgsave(figperf_stab,[expname{1}(2:end),'_laserperfmplotsstab_',savename,'.fig'], '-v7.3');
    saveas(figperf,[expname{1}(2:end),'_laserperfmplotsstab_',savename,'.png']);
    
    hgsave(figperfn_stab,[expname{1}(2:end),'_laserperfmplotsnormstab_',savename,'.fig'], '-v7.3');
    saveas(figperfn_stab,[expname{1}(2:end),'_laserperfmplotsnormstab_',savename,'.png']);
    
    %
    disp('done with saving');
    close all
    clear all
end
