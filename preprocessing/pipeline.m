clear all

%% 
% Large vectors like PdTrig and LaserTrig are not saved. They can easily be
% recalculated later if needed.
%
%
%%
addpath(genpath('/mnt/data/Mitra/cache/codes/matlabkwikreader/'))
addpath(genpath('/mnt/data/Mitra/cache/codes/'))

%% params
make_params_list_VL66
%% constructing paths and assigning channels: be careful that it is /processed or /manual_sorting
% channels
Pdch = n_ephys_channles + n_aux_channels + photodiode_channel;
Laserch = n_ephys_channles + n_aux_channels + laser_channel;

% paths
path2kwe = cell(size(expname));
path2kwd = cell(size(expname));
kwdrecording = ['/recordings/',rec_number,'/data'];
%raw files:
for x=1:1:length(expname)
    path2kwe{x} = ['/mnt/data/Mitra/raw/',animalname,expname{x},'/experiment1.kwe']; % for messages
    path2kwd{x} = ['/mnt/data/Mitra/raw/',animalname,expname{x},'/experiment1_',processor_number,'.raw.kwd']; %for pd % 101 usually
end
% spikes : first experiment in the list only
V1{1}.path2spikes = ['/mnt/data/Mitra/processed/',animalname,expname{1},'_kwd_',processor_number,'_rec',rec_number,'_shkV1_shk0/']; % kilosorted spikes
V1{2}.path2spikes =  ['/mnt/data/Mitra/processed/',animalname,expname{1},'_kwd_',processor_number,'_rec',rec_number,'_shkV1_shk1/'];
LM{1}.path2spikes =  ['/mnt/data/Mitra/processed/',animalname,expname{1},'_kwd_',processor_number,'_rec',rec_number,'_shkLM_shk0/'];
LM{2}.path2spikes =  ['/mnt/data/Mitra/processed/',animalname,expname{1},'_kwd_',processor_number,'_rec',rec_number,'_shkLM_shk1/'];
%%
start_time = cell(size(expname));
rec_length = cell(size(expname));
for x=1:1:length(expname)
    [kwe_infos, messages, ttl] = read_kwe_file(path2kwe{x});
    disp(messages(find(messages.nodeID==100),:).text);
    start_time{x}=messages(find(messages.nodeID==100),:).time_samples
    %messages(find(messages.eventID==0),:).text
    info = h5info(path2kwd{x},kwdrecording);
    rec_length{x}=info.Dataspace.Size(2);
end

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
end_points=cumsum(cell2mat(rec_length));
start_points = [0 end_points(1:end-1)];
%%% MU
LM_Units.MUSpikes = cell(size(expname));
V1_Units.MUSpikes = cell(size(expname));

LM_Units.MUSpikes{1} = LM_Units.MUSpikes_AllExps(find((LM_Units.MUSpikes_AllExps<end_points(1)) )) - start_points(1);
V1_Units.MUSpikes{1} = V1_Units.MUSpikes_AllExps(find((V1_Units.MUSpikes_AllExps<end_points(1)) )) - start_points(1);

for x=2:1:(length(expname))
    LM_Units.MUSpikes{x} = LM_Units.MUSpikes_AllExps(find((LM_Units.MUSpikes_AllExps<end_points(x)) & (LM_Units.MUSpikes_AllExps>end_points(x-1)))) - start_points(x);
    V1_Units.MUSpikes{x} = V1_Units.MUSpikes_AllExps(find((V1_Units.MUSpikes_AllExps<end_points(x)) & (V1_Units.MUSpikes_AllExps>end_points(x-1)))) - start_points(x);
    
end

%%% SU
for i = 1:length(V1_Units.SingleUnitSpikeTimes_AllExps)
    V1_Units.SUSpikes{i} = cell(size(expname));
    V1_Units.SUSpikes{i}{1} = V1_Units.SingleUnitSpikeTimes_AllExps{i}(find((V1_Units.SingleUnitSpikeTimes_AllExps{i}<end_points(1)) )) - start_points(1);
    
    for x=2:1:(length(expname))
        V1_Units.SUSpikes{i}{x} = V1_Units.SingleUnitSpikeTimes_AllExps{i}(find((V1_Units.SingleUnitSpikeTimes_AllExps{i}<end_points(x)) & (V1_Units.SingleUnitSpikeTimes_AllExps{i}>end_points(x-1)))) - start_points(x);
    end
end

for i = 1:length(LM_Units.SingleUnitSpikeTimes_AllExps)
    LM_Units.SUSpikes{i} = cell(size(expname));
    LM_Units.SUSpikes{i}{1} = LM_Units.SingleUnitSpikeTimes_AllExps{i}(find((LM_Units.SingleUnitSpikeTimes_AllExps{i}<end_points(1)) )) - start_points(1);
    
    for x=2:1:(length(expname))
        LM_Units.SUSpikes{i}{x} = LM_Units.SingleUnitSpikeTimes_AllExps{i}(find((LM_Units.SingleUnitSpikeTimes_AllExps{i}<end_points(x)) & (LM_Units.SingleUnitSpikeTimes_AllExps{i}>end_points(x-1)))) - start_points(x);
    end
end
%% parse ssmessages

SSmessages = cell(size(expname));
correcttrial = cell(size(expname));
falsetrial = cell(size(expname));
misstrial = cell(size(expname));
alltrialcf = cell(size(expname));
correcttrialind = cell(size(expname));
falsetrialind = cell(size(expname));
misstrialind = cell(size(expname));
gotrial = cell(size(expname));
nogotrial = cell(size(expname));
alltrialgn = cell(size(expname));
gotrialind = cell(size(expname));
nogotrialind = cell(size(expname));

if do_parse_ssmessages
    ParseSSMessages
end
%% get pd andlaser trigger times
PAllOn = cell(size(expname));
LAllOn = cell(size(expname));
LStepTimeOn = cell(size(expname));
LStepTimeStampOn = cell(size(expname));
LStepTimeStampOff = cell(size(expname));
PStepTimeStampOn = cell(size(expname));
PStepTimeStampOff =cell(size(expname));
LaserTrig = cell(size(expname));
PdTrig = cell(size(expname));
for x=1:1:length(expname)
    
    if exist('SSmessages','var') && numel(SSmessages{x}) > 0
        tStart = min(SSmessages{x}.time_samples) - start_time{x};
        tEnd = max(SSmessages{x}.time_samples) - start_time{x}; 
    else
        tStart = tStart_p;
        tEnd = tEnd_p;
    end
    
    if cut_rec
        tEnd = tEnd_cut;
    end
    
    %%%
    LaserTrig{x} = h5read(path2kwd{x},kwdrecording,[Laserch 1],[1 inf]);
    if strcmp(lasershape ,'ramp')
        [LStepTimeOn{x},LStepTimeStampOn{x}] = FindEdges_rampdown(LaserTrig{x},DiffThL,tStart,tEnd);
    elseif strcmp(lasershape,'square')
       [LStepTimeOn{x},~,LStepTimeStampOn{x},~] = FindEdges(LaserTrig{x},DiffThL,tStart,tEnd);
    else
        disp('choose a lasershape first')
    end
    sprintf('%d laser onset detected',length(LStepTimeStampOn{x}))
    %  sprintf('%d laser offset detected',length(LStepTimeStampOff{x}))
    
    
    
    Window=WindowSec*samplingf;
    LAllOn{x}= repmat(LStepTimeStampOn{x}',[1,2*Window+1])+repmat(double(-Window:Window),[length(LStepTimeStampOn{x}'),1]); % if you don't double, it will be single and the subtraction would be fucked
    pdfig_handle = figure('Units','normalized','Position',[0.2 0.1 0.7 0.7]);
    subplot(3,1,1);plot((-Window:Window)/samplingf,LaserTrig{x}(LAllOn{x}(:,:))');
    %%%
    PdTrig_nofilter = h5read(path2kwd{x},kwdrecording,[Pdch 1],[1 inf]);
    PdTrig{x}=low_pass_filter(PdTrig_nofilter,samplingf,pdfiltfreq);
    [PStepTimeOn,PStepTimeOff,PStepTimeStampOn{x},PStepTimeStampOff{x}] = FindEdges(PdTrig{x},DiffTh,tStart,tEnd);
    
    sprintf('%d pd onset detected',length(PStepTimeStampOn{x}))
    sprintf('%d pd offset detected',length(PStepTimeStampOff{x}))
    
    Window=WindowSec*samplingf;
    WindowL=WindowSec*samplingf;
    PAllOn{x}= repmat(PStepTimeStampOn{x}',[1,(Window+WindowL)+1])+repmat(double(-WindowL:Window),[length(PStepTimeStampOn{x}'),1]); % if you don't double, it will be single and the subtraction would be fucked
    subplot(3,1,2);plot((-WindowL:Window)/samplingf,PdTrig{x}(PAllOn{x}(:,:))');
    subplot(3,1,3);plot((-WindowL:Window)/samplingf,LaserTrig{x}(PAllOn{x}(:,:))');
    pdfig_handle.Children(3).Title.String = sprintf('exp number %d\n%d pd onset detected\n%d pd offset detected\n%d laser onset detected',...
        x,length(PStepTimeStampOn{x}),length(PStepTimeStampOff{x}),length(LStepTimeStampOn{x}));
    pdfig_handle.Children(3).Title.FontWeight = 'normal';
    pdfig_handle.Children(3).Title.FontSize = 9;

    figure;plot((-WindowL:Window)/samplingf,PdTrig_nofilter(PAllOn{x}(:,:))')
end
%% behavior analysis

if do_behavior_analysis
    %simple_behavior;
    speed_revised;
    isrevisedspeed=1;
    
    fprintf('Rslag = %d\n', cell2mat(Rslag))
    x=1;
    fprintf('number of no grooming trials= %d\n', numel(nogroomingind{x}))   
end
%% binning laser delays
x=1;
[laserdelayfig_handle,LaserDelay,centers,LaserDelayBinned] = Extractlaserdelaynan(LStepTimeOn{x},PAllOn{x},Window,samplingf,Laserdelaybinwindow,Laserdelaybinnum);
%% separate v1 units that are excited by laser - done by binning
%% saving
prompt = 'Are you happy with preprocessing results? Y/N [Y]: ';
userin = input(prompt,'s');
if strcmp(userin,'y') || strcmp (userin,'Y')
    cd(savingdirroot)
    x=1; % now, only saving the first experiment
    
    % make folders if necessary
    if ~exist(fullfile(savingdirroot,'preprocessing'),'dir')
        mkdir('preprocessing')
    end
    % timestamped name
    cd('preprocessing')
    savename=datestr(now,'yyyy_mm_dd_HH_MM_SS');
    mkdir(savename);
    cd(savename);
    % save variables
    save([expname{x}(2:end),'_preprocessing_',savename],...
        'animalname','V1','V1_Units','LM','LM_Units','units_to_extract',...
        'Window','WindowL','WindowSec','ttl','tStart','tEnd','tStart_p','tEnd_p','start_time','start_points',...
        'savingdirroot','samplingf','rec_number','rec_length','path2kwd','path2kwe','expname','DiffTh','DiffThL',...
        'processor_number','photodiode_channel','Netevent_nodeID','n_ephys_channles','n_aux_channels','Pdch',...
        'Laserch','laser_channel','kwdrecording','kwe_infos','end_points','do_behavior_analysis',...
        'do_parse_ssmessages','messages','LAllOn','PAllOn','pdfiltfreq','Laserdelaybinnum','Laserdelaybinwindow',...
        'LStepTimeOn','LStepTimeStampOff','LStepTimeStampOn','centers','LaserDelay','LaserDelayBinned',...
        'PStepTimeOff','PStepTimeOn','PStepTimeStampOff','PStepTimeStampOn',...
        'SSmessages','correcttrial','falsetrial','misstrial','alltrialcf','correcttrialind',...
        'falsetrialind','misstrialind','gotrial','nogotrial','alltrialgn','gotrialind','nogotrialind',...
        'lick_channel', 'running_channel','grooming_speed_threshold','grooming_lick_threshold',...
        'isrevisedspeed','speedWindowms', 'speedticks', 'speed','licks', 'nogroomingind', 'cut_rec', 'tEnd_cut','-v7.3');
    
    %save figures
   % saveas(pdfig_handle,[expname{x}(2:end),'_pd_',savename,'.fig']);
    saveas(pdfig_handle,[expname{x}(2:end),'_pd_',savename,'.png']);
    
   % saveas(laserdelayfig_handle,[expname{x}(2:end),'_laserdelay_',savename,'.fig']);
    saveas(laserdelayfig_handle,[expname{x}(2:end),'_laserdelay_',savename,'.png']);
    
    if do_behavior_analysis
        % for some reason, saving figs errors
        
        %saveas(beh_f1,[expname{x}(2:end),'_behavior_',savename,'.fig']);
        saveas(beh_f1,[expname{x}(2:end),'_behavior_',savename,'.png']);
    
        %saveas(beh_f2,[expname{x}(2:end),'_behaviorexgrooming_',savename,'.fig']);
        saveas(beh_f2,[expname{x}(2:end),'_behaviorexgrooming_',savename,'.png']);
   
    end
    
    %
    disp('done with saving');
    close all
    clear all   
end
