% per animal, we make stimulus induces lfp and csd (from saved downsampe lfp code)
% csd will be used later to map channel number to depth
% to be saved as a separate .mat file per animal
% if spatially downsampled, lfp is 1 per 4 channels. so 8 channels per shank (~100
% micrometer distance between)
% output:
%      - data_lfp_bs, struct with go/ nogo fields, each a tensor of
%       nchannels*ntrials* time. LFP for baseline trials
%      - data_lfp: same as above but for different laser lags
%
% Before running, make sure the animal has the lfp file, task file, and
% shankinfo file:
% order = {'V1_0','V1_1','LM_0','LM_1'};
% save('shankinfo.mat','order');
%
% !!! IMPORTANT: depth is based on csd
%%
clear all
%%
addpath(genpath('/mnt/data/Mitra/cache/codes/matlabkwikreader/'))
addpath(genpath('/mnt/data/Mitra/cache/codes/'))
%%
dsrate = 50;
nchannels = 128;    %for all 4 shanks (can't do only 1 shank, order file gets corrupted  -
% and the same value for all animals (32 for 1:4 spatially downsampled lfp, 128 full)
savevars = 1;
lfptimewindow = 30*[-100 600];%
smoothinglength = 5;
%% start making stimulus aligned LFP
 
%  animallist = {'MPV17','MPV18_2',...
%      'VL53','VL52','VL51','VL66'};%,...
% %        'VL61','VL63','VL55','VL59','VL50'};
%  preprocessinglist = {'2019_04_05_21_25_20','2019_04_05_15_04_39',...
%      '2017_12_05_18_58_17','2017_12_04_15_08_47','2017_12_03_17_30_15','2018_03_20_11_53_19',...
%       '2018_03_16_16_53_13','2018_04_01_12_24_38','2018_04_01_15_03_48','2018_04_01_13_24_48','2017_12_01_16_21_49'};

animallist = {'MPV32_2','MPV33','MPV31','MPV34_2',...
    'MPV36_2','MPV30','MPV35_2'};
preprocessinglist = {'2020_01_21_18_13_31','2020_01_21_17_12_32','2020_01_22_13_47_56','2020_01_22_14_28_52',...
    '2020_01_22_11_38_57','2020_01_22_14_51_15','2020_01_28_20_29_59'};



trialtypelist = {'gotrialind','nogotrialind'};
targettrialtypelist = {'correctgotrialind','correctnogotrialind'};
trialtype_namelist = {'go','nogo'};


for animal_i=1:length(animallist) % 1
    cd('/mnt/data/Mitra/figs');
    cd(sprintf('%s/preprocessing/%s',animallist{animal_i},preprocessinglist{animal_i}))
    
    matfilename = dir('task_*.mat');
    lfpfilename = dir('onlytime_downsampled*.mat');
    shankfilename = dir('shanki*.mat');
    load(matfilename.name);
    load(lfpfilename.name);
    load(shankfilename.name);
    clear valves licks allvalve alllick PStepTimeOn PStepTimeOff LStepTimeOn
    %%%%%%%%%%%%%%%%%% make shank order table:
    animalshankorder = cell(1,nchannels);
    for i=1:nchannels
        if i<=nchannels/4
            animalshankorder{i} = order{1};
        elseif i<=nchannels/2
            animalshankorder{i} = order{2};
        elseif i<=3*nchannels/4
            animalshankorder{i} = order{3};
        else
            animalshankorder{i} = order{4};
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%% define params after loading mat files
    clear stabletrialind
    min_n_trials_per_delay = 10;%10 when excludegrooming=0. set to 15 for FF
    %
    excludegrooming = 1; %0 -1:only grooming. When this is -1, only correct should be 0. For correct trials, it is 1
    onlycorrect = 1; % change to 1 for correct trials only.
    onlystableperiod = 1;
    %%%%%%%
    
    timepoints = 1:1:rec_length;
    timepoints_ds = 1:dsrate:rec_length;
    
    data_lfp_bs = struct;
    data_lfp = struct;
    
    for trialtype_i = 1:2
        trialtype = trialtypelist{trialtype_i};
        targettrialtype = targettrialtypelist{trialtype_i};
        trialtype_name = trialtype_namelist{trialtype_i};
        %%%%%%%%% basline
        % indices of baseline trials
        if excludegrooming == 1
            indss0=intersect(intersect(find(isnan(LaserDelay)),eval(trialtype)),nogroomingind);
        elseif excludegrooming == 0
            indss0= intersect(find(isnan(LaserDelay)),eval(trialtype));
        elseif excludegrooming == -1
            indss0=setdiff(intersect(find(isnan(LaserDelay)),eval(trialtype)),nogroomingind);
        end
        if onlycorrect
            indss0=intersect(indss0,eval(targettrialtype));
        end
        if onlystableperiod & exist('stabletrialind','var')
            indss0=intersect(indss0,stabletrialind);
        end
        
        PAllOn_bs = PAllOn(indss0,(floor(size(PAllOn,2)/2)+lfptimewindow(1)):(floor(size(PAllOn,2)/2)+lfptimewindow(2)));
        data_lfp_bs_temp =single(nan([nchannels, size(PAllOn_bs)]));
        for ch=1:1:nchannels
            for indss0_i = 1:length(indss0)
                querytimepoints = double(PAllOn_bs(indss0_i,:));
                data_lfp_bs_temp(ch,indss0_i,:) = single(interp1(timepoints_ds,chdata_ds(ch,:),querytimepoints));
            end
        end
        %%%%%%%%%%% laserdelay
        data_lfp_temp = cell(0,0);
        for ldelay =1:length(centers)
            ldelay
            indss= intersect(find(LaserDelayBinned == (ldelay)),eval(trialtype));
            if length(indss) > min_n_trials_per_delay
                if excludegrooming == 1
                    indss=intersect(indss,nogroomingind);
                elseif  excludegrooming == -1
                    indss=setdiff(indss,nogroomingind);
                end
                if onlycorrect
                    indss=intersect(indss,eval(targettrialtype));
                end
                if onlystableperiod & exist('stabletrialind','var')
                    indss=intersect(indss,stabletrialind);
                end
                
                
                PAllOn_l = PAllOn(indss,(floor(size(PAllOn,2)/2)+lfptimewindow(1)):(floor(size(PAllOn,2)/2)+lfptimewindow(2)));
                data_lfp_temp{end+1} =single(nan([nchannels, size(PAllOn_l)]));
                for ch=1:1:nchannels
                    for indss_i = 1:length(indss)
                        querytimepoints = double(PAllOn_l(indss_i,:));
                        data_lfp_temp{end}(ch,indss_i,:) = single(interp1(timepoints_ds,chdata_ds(ch,:),querytimepoints));
                    end
                end
            end
            
        end
        %%%%%%%%%%%%%%%%%%%%
        data_lfp_bs = setfield(data_lfp_bs,trialtype_name,data_lfp_bs_temp);
        data_lfp = setfield(data_lfp,trialtype_name,data_lfp_temp);
    end
    % free up memory
    clear chdata_ds  data_lfp_bs_temp data_lfp_temp
    %% get csd and lfp for go/nogo baseline - get depth from L4 border
    % function to calculate depth based on go stimulus aligned lfp

    % the last (deepest) channel at the end of the high blue blob - still blue (for V1
    % is L4, for LM, would just be  a ref)
    [shankdepthmap,refdepth,lfpdepth_f] = alignlfpandgetdepth(order,data_lfp_bs,animalshankorder,lfptimewindow,smoothinglength);
    %% saving
    if savevars
        cd('/mnt/data/Mitra/figs');
        cd(sprintf('%s/preprocessing/%s',animallist{animal_i},preprocessinglist{animal_i}))
        % save depth to shankinfo.mat
        save('shankinfo.mat','refdepth','shankdepthmap','-append');
        % save stimulus aligned lfp:
        save('stimalignedlfp.mat','data_lfp_bs','data_lfp','lfptimewindow','-v7.3');
        % save lfpdepth figures
        for shanknum = 1:numel(order)
            saveas(lfpdepth_f{shanknum},sprintf('lfp_csd_depth_%s.fig',order{shanknum}))
        end
        close all
    end
end
%% example plots of go/nogo lfp in baseline vs laser trials :
figure;
im1= imagesc([squeeze(nanmean(data_lfp_bs.go,2)),squeeze(nanmean(data_lfp_bs.nogo,2))]);
im1.Parent.CLim = [-2000 2000];
figure;
im1= imagesc([squeeze(nanmean(data_lfp.go{2},2)),squeeze(nanmean(data_lfp.nogo{2},2))]);
im1.Parent.CLim = [-2000 2000];
