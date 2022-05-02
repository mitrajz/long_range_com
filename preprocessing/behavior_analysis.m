


animallist = {'MPV14_2','MPV14_2','MPV17','MPV17',...
    'MPV18','MPV18_2','MPV16_2'};
preprocessinglist = {'2019_03_08_15_20_04','2019_03_14_12_02_04','2019_03_03_13_18_55','2019_03_12_13_34_38',...
    '2019_02_28_21_05_57','2019_03_03_15_58_19','2019_03_08_13_35_48'};


for animal_i=1:length(animallist)
    
    x=1;
    speedWindowms =10;
    checkpos = 0;
    isrevisedspeed = 0;
    speedticks = nan;
    speed =nan;
    
    animallist{animal_i}
    cd('/mnt/data/Mitra/figs');
    cd(sprintf('%s/preprocessing/%s',animallist{animal_i},preprocessinglist{animal_i}))
    matfilename = dir('task*.mat');
    load(matfilename.name)
    
    
    valvech = n_ephys_channles + n_aux_channels + 5;
    lickch = n_ephys_channles + n_aux_channels + 6;
    RSch=n_ephys_channles + n_aux_channels + running_channel;
    
    if checkpos
        Pos = h5read(path2kwd{x},kwdrecording,[RSch 1],[1 inf]);
        Rslag{x} = find(Pos > 0.99*max(Pos),1) - find(Pos < 1.1*min(Pos),1);
        
        if Rslag{x}<15000 || Rslag{x}>45000
            Rslag{x} = 30000;
        end
    end
    
    gain=h5read(path2kwd{x},'/recordings/0/application_data/channel_bit_volts');
    
    %%% lick and valve
    alllick{x} = h5read(path2kwd{x},kwdrecording,[lickch 1],[1 inf]);
    allvalve{x} = h5read(path2kwd{x},kwdrecording,[valvech 1],[1 inf]);
    
    %% normalize lick and valves
    lickthreshold = prctile(alllick{x},95);
    valvethreshold = prctile(allvalve{x},99.5);
    
    lick_th = alllick{x};
    lick_th(alllick{x}<lickthreshold) = 0;
    lick_th(alllick{x}>=lickthreshold) = 1;
    
    valve_th = allvalve{x};
    valve_th(allvalve{x}<valvethreshold) = 0;
    valve_th(allvalve{x}>=valvethreshold) = 1;
    
    %%%
    if checkpos
        speedticks = nan(size(PAllOn{x},1),size(PAllOn{x},2)-1);
        numspeedbins = (Window+WindowL)/(speedWindowms*(samplingf/1000));
        speed= nan(size(PAllOn{x},1),numspeedbins);
    end
    licks = nan(size(PAllOn{x}));
    valves = nan(size(PAllOn{x}));
    
    for trialind=1:1:size(PAllOn{x},1)

        if checkpos
            if ((PAllOn{x}(trialind,end) + Rslag{x}) - length(Pos)) <= 0
                %%%% Rs
                trialPos=double(Pos(PAllOn{x}(trialind,:) + Rslag{x})')*double(gain(RSch));
                %figure;plot(trialPos)
                speedticks(trialind,:)=double(diff(trialPos));
                speedticks(trialind,(speedticks(trialind,:))<0)=nan; % negative: going back or resetting position
                speedticks(trialind,(speedticks(trialind,:)<0.01))=0; % noise around a contant value
                speedticks(trialind,(speedticks(trialind,:)>0.5)) =nan; % >0.5cm steps
                
                for n=1:numspeedbins
                    speed(trialind,n) = nansum(speedticks(trialind,(1+(n-1)*(speedWindowms*(samplingf/1000))):(n*(speedWindowms*(samplingf/1000)))))...
                        /(speedWindowms/1000);% in cm/s
                end
                %figure;plot(speed)
                %figure;plot(RSs{x}(PAllOn{x}(trialind,:) + Rslag{x})')
                %%%% licks
            end
        end
        licks(trialind,:) = lick_th(PAllOn{x}(trialind,:));
        valves(trialind,:) = valve_th(PAllOn{x}(trialind,:));
        
    end
    
    %% first licks and comparison to valve
    lickcontforms = 10;
    lickcontforsample = lickcontforms*30;
    firstlicksample = nan(1,size(PAllOn{x},1));
    firstvalvesample = nan(1,size(PAllOn{x},1));
    thr = 0.2; 
    for trialind=1:1:size(PAllOn{x},1)

        fl=find(single(smooth(licks(trialind,floor(size(licks,2)/2):end),lickcontforsample)) > thr,1)-lickcontforsample/2;
        fv=find(single(valves(trialind,floor(size(licks,2)/2):end)) == 1,1);
        if numel(fl) && fl>0
            firstlicksample(trialind) = fl;
        end
        
        if numel(fv) && fv>0
            firstvalvesample(trialind) = fv;
        end
    end
    f_lick1=figure;plot(firstlicksample,'.')
    hold on;plot(firstvalvesample,'r.')
    f_lick2=figure;histogram(firstlicksample/30, 100,'Normalization','count','FaceColor','b')
    hold on;histogram(firstvalvesample/30, 100,'Normalization','count','FaceColor','r')
  
    respwindowms = [80 prctile(firstvalvesample,99)/30];
    repwindowsample = 30*respwindowms;
    

    
    correctnogotrialind{1} = [];
    correctgotrialind{1} = [];
    incorrectnogotrialind{1} = [];
    incorrectgotrialind{1} = [];
    
    FA=0;
    Miss=0;
    
    indss_nogo = nogotrialind{x};
    indss_go = gotrialind{x};

    for i = 1:length(indss_nogo)
        if (firstlicksample(indss_nogo(i))>repwindowsample(1)) & (firstlicksample(indss_nogo(i))<repwindowsample(2))
            FA=FA+1;
            incorrectnogotrialind{1}(end+1) = indss_nogo(i);
        else
            correctnogotrialind{1}(end+1) = indss_nogo(i);
        end
    end
    
    for i = 1:length(indss_go)
        if (firstlicksample(indss_go(i))>repwindowsample(1)) & (firstlicksample(indss_go(i))<repwindowsample(2))
            correctgotrialind{1}(end+1) = indss_go(i);
        else
            Miss=Miss+1;
            incorrectgotrialind{1}(end+1) = indss_go(i);
        end
    end
    
    FA= 100*FA/length(indss_nogo)
    Miss=100*Miss/length(indss_go)
    
    correcttrials = zeros(1,length(size(PAllOn{x},1)));
    correcttrials([correctnogotrialind{1} correctgotrialind{1}]) = 1;
    f_behovertime = figure;plot(movmean(correcttrials, [0 50]));
    %% grooming

    nogroomingind=cell(size(expname));
    nogroomingind{x} = [];
    for trialind=1:1:size(PAllOn{x},1) % grooming calculated only for prestimulus window
        if ((nanmean(licks(trialind,1:floor(size(licks,2)/2))) < grooming_lick_threshold) )
            nogroomingind{x}(end+1) = trialind;
        end
    end

    %% plots
    
    figure;
    beh_f1=gcf;beh_f1.Name= 'including grooming';
    imagesc(((-Window+1)/samplingf:(speedWindowms/1000):Window/samplingf),[],...
        [licks(gotrialind{x},:);nan(1,size(licks,2));licks(nogotrialind{x},:)]);
    hold on;
    ax1=imagesc(((-Window+1)/samplingf:(speedWindowms/1000):Window/samplingf),[],...
        [valves(gotrialind{x},:);nan(1,size(valves,2));valves(nogotrialind{x},:)],'AlphaData',0.5);
    colormap(ax1.Parent,'gray')
    ax1.Parent.XLabel.String='time(s)';
    ax1.Parent.YLabel.String='trials';
    ax1.Parent.Title.String='licks aligned to stimulus onset';
    
    
    figure;
    beh_f2=gcf;beh_f2.Name= 'excluding grooming';
    imagesc(((-Window+1)/samplingf:(speedWindowms/1000):Window/samplingf),[],...
        [licks(intersect(gotrialind{x},nogroomingind{x}),:);nan(1,size(licks,2));licks(intersect(nogotrialind{x},nogroomingind{x}),:)]);
    hold on;
    ax1=imagesc(((-Window+1)/samplingf:(speedWindowms/1000):Window/samplingf),[],...
        [valves(intersect(gotrialind{x},nogroomingind{x}),:);nan(1,size(valves,2));valves(intersect(nogotrialind{x},nogroomingind{x}),:)],'AlphaData',0.5);
    colormap(ax1.Parent,'gray')
    ax1.Parent.XLabel.String='time(s)';
    ax1.Parent.YLabel.String='trials';
    ax1.Parent.Title.String='licks aligned to stimulus onset';
    %% saving
    x=1; % only saving the first experiment
    % timestamped name
    savename=datestr(now,'yyyy_mm_dd_HH_MM_SS');
    % save variables
    % beh file
    save(['beh_',expname{x}(2:end),'_preprocessing_',savename],...
        'animalname','grooming_lick_threshold','nogroomingind',...
        'correcttrials','FA','Miss','correctnogotrialind','incorrectnogotrialind','correctgotrialind','incorrectgotrialind',...
        'isrevisedspeed','speedWindowms', 'speedticks', 'speed',...   
        'firstlicksample','firstvalvesample','lickcontforms','thr','valves','licks','lickthreshold','valvethreshold',...
        'allvalve','alllick','checkpos','-v7.3');
    % preprocesing file
    % remove this for animals with later preprocessing
    save([expname{x}(2:end),'_preprocessing_',savename],...
        'animalname','V1','V1_Units','LM','LM_Units','units_to_extract',...
        'Window','WindowL','WindowSec','ttl','tStart','tEnd','tStart_p','tEnd_p','start_time','start_points',...
        'savingdirroot','samplingf','rec_number','rec_length','path2kwd','path2kwe','expname','DiffTh','DiffThL',...
        'processor_number','photodiode_channel','Netevent_nodeID','n_ephys_channles','n_aux_channels','Pdch',...
        'Laserch','laser_channel','kwdrecording','kwe_infos','end_points','do_behavior_analysis',...
        'do_parse_ssmessages','messages','LAllOn','PAllOn','pdfiltfreq','Laserdelaybinnum','Laserdelaybinwindow',...
        'LStepTimeOn','LStepTimeStampOff','LStepTimeStampOn','centers','LaserDelay','LaserDelayBinned',...
        'PStepTimeOff','PStepTimeOn','PStepTimeStampOff','PStepTimeStampOn',...
        'SSmessages',...
        'gotrial','nogotrial','alltrialgn','gotrialind','nogotrialind',...
        'lick_channel', 'running_channel','grooming_lick_threshold',...
        'cut_rec', 'tEnd_cut','-v7.3');
    
    %save figures
    
    saveas(beh_f1,[expname{x}(2:end),'_beh1_',savename,'.png']);
    saveas(beh_f2,[expname{x}(2:end),'_beh2_',savename,'.png']);
    saveas(f_lick1,[expname{x}(2:end),'_lick1_',savename,'.png']);
    saveas(f_lick2,[expname{x}(2:end),'_lick2_',savename,'.png']);
    saveas(f_behovertime,[expname{x}(2:end),'_perfmtime_',savename,'.png']);
    
    clearvars -except animal_i animallist preprocessinglist
    close all
end