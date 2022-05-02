%%
clear all
%%
addpath(genpath('/mnt/data/Mitra/cache/codes/matlabkwikreader/'))
addpath(genpath('/mnt/data/Mitra/cache/codes/'))

% animallist = {'MPV17','MPV18_2',...
%     'VL53','VL52','VL51','VL66',...
%     'VL61','VL63','VL55','VL59','VL50'};
% preprocessinglist = {'2019_04_05_21_25_20','2019_04_05_15_04_39',...
%     '2017_12_05_18_58_17','2017_12_04_15_08_47','2017_12_03_17_30_15','2018_03_20_11_53_19',...
%     '2018_03_16_16_53_13','2018_04_01_12_24_38','2018_04_01_15_03_48','2018_04_01_13_24_48','2017_12_01_16_21_49'};

animallist = {'MPV32_2','MPV33','MPV31','MPV34_2',...
    'MPV36_2','MPV30','MPV35_2'};
preprocessinglist = {'2020_01_21_18_13_31','2020_01_21_17_12_32','2020_01_22_13_47_56','2020_01_22_14_28_52',...
    '2020_01_22_11_38_57','2020_01_22_14_51_15','2020_01_28_20_29_59'};

pathchangenecessary = 0;

for animal_i=1:length(animallist)
    cd('/mnt/data/Mitra/figs');
    cd(sprintf('%s/preprocessing/%s',animallist{animal_i},preprocessinglist{animal_i}))
    
    matfilename = dir('task_*.mat');
    load(matfilename.name);
    clearvars -except path2kwd kwdrecording rec_length animal_i animallist preprocessinglist pathchangenecessary
    
    if pathchangenecessary
        if animal_i < 3
            %MPV:
            % reading for network
            % path2kwd=['/mnt/javadzam/winstor/swc/mrsic_flogel/public',path2kwd(16:end)];
            % or locally
            path2kwd=['/mnt/data/Mitra/raw',path2kwd(77:end)];
        else
            % %VL
            % path2kwd{1}=['/mnt/javadzam/winstor/swc/mrsic_flogel/public/projects/MiJa_20160601_VisualLongRangeConnectivity/Ephys',...
            %     path2kwd{1}(16:end)];
            path2kwd=['/mnt/data/Mitra/raw',path2kwd(20:end)];
        end
    end
    
    %% params
    downsampligmode = 2; % 1: downsample in both time and space 2: downsample only in time
    nchannels=128;
    dsrate = 50;% 30000 to 600
    save_result = 1;
    %%
    if downsampligmode == 1
        sp_dsrate = 4;
    elseif downsampligmode == 2
        sp_dsrate = 1;
    end
   
    %chdata = zeros(1,rec_length);
    chdata_ds = zeros(nchannels,ceil(rec_length/dsrate));
    %%
    %%% pre filter: below 300
    [b,a] = butter(5,(300/(30000/2)),'low');
    for i=1:sp_dsrate:nchannels
        animallist{animal_i}
        i
        tic
        chdata = double(h5read(path2kwd,kwdrecording,[i 1],[1 inf]));
        chdata=filtfilt(b,a, chdata);
        %     %%% and notch 50
        %     [b,a] = butter(3,([47 53]/(30000/2)),'stop');
        %     chdata=filtfilt(b,a, chdata);
        chdata_ds(i,:) = chdata(1:dsrate:end);
        toc
    end
    chdata_ds(find(mean(chdata_ds,2)==0),:)=[];
    %%
    if save_result
        cd('/mnt/data/Mitra/figs');
        cd(sprintf('%s/preprocessing/%s',animallist{animal_i},preprocessinglist{animal_i}))
        save('onlytime_downsampled_lfp.mat','chdata_ds','dsrate','-v7.3')
    end
end
