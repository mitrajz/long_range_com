% make generic LM and V1 cells (both LM and V1), static mice (running speed = 0)
clear all
%%
addpath(genpath('/mnt/data/Mitra/cache/codes/matlabkwikreader/'))
addpath(genpath('/mnt/data/Mitra/cache/codes/'))
%%

LMcells = cell(0,0);
V1cells = cell(0,0);
Behcells = cell(0,0);


animallist = {'MPV17','MPV18_2',...
    'VL53','VL52','VL51','VL66'};%,...
    %'VL61','VL63','VL55','VL59','VL50'};
preprocessinglist = {'2019_04_05_21_25_20','2019_04_05_15_04_39',...
    '2017_12_05_18_58_17','2017_12_04_15_08_47','2017_12_03_17_30_15','2018_03_20_11_53_19',...
    '2018_03_16_16_53_13','2018_04_01_12_24_38','2018_04_01_15_03_48','2018_04_01_13_24_48','2017_12_01_16_21_49'};
exptype = {'FB','FB',...
     'FB','FB','FB','FB',...
     'FF','FF','FF','FF','FF'};


% rotate


% animallist = {'MPV17','MPV18_2'};
% preprocessinglist = {'2019_04_07_16_57_46','2019_04_05_16_48_40'};
%  exptype = {'FB','FB','FB'};
% animallist = {'MPV18_2'};
% preprocessinglist = {'2019_04_05_16_48_40'};
% exptype = {'FB','FB','FB'};


% quick performance check: doesn't work when AR was 1. Should calculate
% first lick
%[size(intersect(misstrialind,nogroomingind))/size(intersect(gotrialind,nogroomingind)) size(intersect(falsetrialind,nogroomingind))/size(intersect(nogotrialind,nogroomingind))  ]

%
trialtypelist = {'gotrialind','nogotrialind'};
targettrialtypelist = {'correctgotrialind','correctnogotrialind'};
trialtype_namelist = {'go','nogo'};


for animal_i=1:length(animallist) % 1
    cd('/mnt/data/Mitra/figs');
    cd(sprintf('%s/preprocessing/%s',animallist{animal_i},preprocessinglist{animal_i}))
    if animal_i<3
        matfilename = dir('task_*.mat');
    else
        matfilename = dir('s*.mat');
        perfmatfilename = dir('p*.mat');
        load(perfmatfilename.name)
    end
    load(matfilename.name)
    %%%%%%%%%%%%%%%%%%%%%%%%define params after loading mat files
    clear stabletrialind
    edgestep = 80*30;
    min_n_trials_per_delay = 10;%10 when excludegrooming=0. set to 15 for FF
    nbins = 1;%1
    %
    excludegrooming = 1; %0 -1:only grooming. When this is -1, only correct should be 0. For correct trials, it is 1
    onlycorrect = 1; % change to 1 for correct trials only.
    onlystableperiod = 1;
    %%%%%%%%%%%%%%%%%%%%%%%%
    

    
    if animal_i>2
        PAllOn = PAllOn{1};
        %licks = licks{1};
        nogroomingind = nogroomingind{1};
        nogotrialind = nogotrialind{1};
        gotrialind = gotrialind{1};      
        correctgotrialind = correctgotrialind{1};
        correctnogotrialind = correctnogotrialind{1};
        incorrectgotrialind = incorrectgotrialind{1};
        incorrectnogotrialind = incorrectnogotrialind{1};
    end
    
    for area = {'V1','LM'}
        for shank=1:2
            for neuronnum =1:size(eval(sprintf('%s{shank}.SingleUnitSpikeTimes',area{:})),2)
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%% make thiscell struct
                clear thiscell
                % trialtype independent
                thiscell.animal = animalname;
                thiscell.simulcode = animal_i;
                thiscell.exptype = exptype{animal_i};
                thiscell.shank = shank;
                thiscell.neuronnum = neuronnum;
                thiscell.expname = expname;
                thiscell.edgestep = edgestep;
                thiscell.depth = eval(sprintf('%s{shank}.SingleUnitsDepths(neuronnum)',area{:}));
                %      thiscell.spiketemplate = eval(sprintf('%s{shank}.SingleUnitTemplates{neuronnum}',area{:}));
                if animal_i<3
                    thiscell.drift = eval(sprintf('%s{shank}.SingleUnitDriftMeasure_stable(neuronnum)',area{:}));
                end
                %%%%%%% trialtype dependent
                % empty structs
                thiscell.nbs = struct;
                thiscell.nls = struct;
                thiscell.nlsAv = struct;
                thiscell.nbsAv = struct;
                thiscell.smb = struct;
                thiscell.smb_p = struct;
                thiscell.smb_z = struct;
                thiscell.smb_bs = struct;
                thiscell.smb_ls = struct;
                thiscell.smb_delta = struct;
                thiscell.smb_centers = struct;
                thiscell.averagesilencing = struct;
                thiscell.laAbs = struct;
                thiscell.laAls = struct;
                %%%%%%%%%%%%%%%%%%%%%%%%%
                thisBeh.licks = struct;
                thisBeh.licksalltrials = struct;
                thisBeh.ldelayindss = struct;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                
                for trialtype_i = 1:2
                    trialtype = trialtypelist{trialtype_i};
                    targettrialtype = targettrialtypelist{trialtype_i};
                    trialtype_name = trialtype_namelist{trialtype_i};
                    %%%%%%%%%
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
                    %%%%
                    %%%%%% response of the cell in baseline trials
                    
                    units=eval(sprintf('%s{shank}.SingleUnitSpikeTimes{neuronnum}',area{:}));
                    rasterbinnedn=nan(size(PAllOn,1),(length(1:edgestep:size(PAllOn,2))-1));
                    rasterbinnedlicks=nan(size(PAllOn,1),(length(1:edgestep:size(PAllOn,2))-1));
                    for i=1:1:size(PAllOn,1)
                        edges = PAllOn(i,1):edgestep:PAllOn(i,end);
                        [rasterbinnedn(i,:),~] = histcounts(units,edges);
                        edgeslick = 0:edgestep:size(licks,2)-1;
                        [rasterbinnedlicks(i,:),~] = histcounts(find(licks(i,:)),edgeslick);
                    end
                    Avfr = samplingf*length(units)/double(units(end) - units(1));
                    Depth = eval(sprintf('%s{shank}.SingleUnitsDepths(neuronnum)',area{:}));
                    %%%%%% silencing effect
                    clear smb stemb
                    smb = nan(1,max(unique(LaserDelayBinned)));
                    smb_p = nan(1,max(unique(LaserDelayBinned)));
                    smb_z = nan(1,max(unique(LaserDelayBinned)));
                    smb_bs = nan(1,max(unique(LaserDelayBinned)));
                    smb_ls = nan(1,max(unique(LaserDelayBinned)));
                    stemb = nan(1,max(unique(LaserDelayBinned)));
                    smb_centers = nan(1,max(unique(LaserDelayBinned)));
                    laAbs = cell(0,0);
                    laAls = cell(0,0);
                    nls = cell(0,0);
                    
                    smb_delta = nan(1,max(unique(LaserDelayBinned)));
                    stemb_delta = nan(1,max(unique(LaserDelayBinned)));
                    
                    nlsAv = nan(length(centers),size(rasterbinnedn,2));
                    
                    ldelayindss = cell(0,0);
                    ldelayindss{end+1} = indss0;
                    for ldelay =1:length(centers)

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
                            
                            ldelayindss{end+1} = indss;
                            nlsAv(ldelay,:) = nanmean(rasterbinnedn(indss,:),1);
                            
                            %%%%%%%%%%%%%%%%%% laser aligned
                            rasterbinned_la=nan(size(PAllOn,1),1);
                            for i=1:1:size(PAllOn,1)
                                point=floor(floor(size(PAllOn,2)/2)+centers(ldelay)*30);
                                edges =PAllOn(i,point):edgestep:PAllOn(i,point+1*edgestep);
                                [rasterbinned_la(i,:),~] = histcounts(units,edges);
                            end
                            
                            laAbs{end+1} = rasterbinned_la(indss0);
                            laAls{end+1} = rasterbinned_la(indss);
                            nls{end+1} = rasterbinnedn(indss,:);
                            %%%%%%%%%%%%%%%%%%
                            
                            
                            if false
                                
                                
                                %convert centers(ldelay) to bins
                                duration = (size(rasterbinnedn,2)/2+1+ceil(centers(ldelay)/(edgestep/30)):size(rasterbinnedn,2)/2+1+ceil(centers(ldelay)/(edgestep/30))+nbins);
                                
                                
                                if mean(sum(rasterbinnedn(indss0,duration),2)) > 0.1 %
                                    differencefromavbaseline = 100*(sum(rasterbinnedn(indss,duration),2) - mean(sum(rasterbinnedn(indss0,duration),2))) / mean(sum(rasterbinnedn(indss0,duration),2));
                                    deltaf = (sum(rasterbinnedn(indss,duration),2) - mean(sum(rasterbinnedn(indss0,duration),2)));
                                    
                                    smb_delta(ldelay) = mean(deltaf,1);
                                    if numel(sum(rasterbinnedn(indss,duration),2))
                                        [smb_p(ldelay),~] =  ranksum(sum(rasterbinnedn(indss0,duration),2),sum(rasterbinnedn(indss,duration),2));
                                    end
                                    smb_z(ldelay) = (nanmean(sum(rasterbinnedn(indss0,duration),2)) - nanmean(sum(rasterbinnedn(indss,duration),2)))/...
                                        sqrt((nanmean(sum(rasterbinnedn(indss0,duration),2))/length(indss0))+(nanmean(sum(rasterbinnedn(indss,duration),2))/length(indss)));
                                    
                                    smb(ldelay) = mean(differencefromavbaseline,1) ;
                                    smb_bs(ldelay) = mean(mean(sum(rasterbinnedn(indss0,duration),2)),1);
                                    %stdmb = std(differencefromavbaseline,0,1);
                                    stemb(ldelay) = std(differencefromavbaseline,1)/sqrt(length(indss));
                                    stemb_delta(ldelay) = std(deltaf,1)/sqrt(length(indss));
                                    
                                    
                                    
                                end
                            end
                            %%%%%%%%%%%%%%%%%%%%%%%%% with laser aligned
                            
                            %convert centers(ldelay) to bins
                            duration = 1;
                         
                            if nanmean(sum(rasterbinned_la(indss0,duration),2)) > 0.2
                                differencefromavbaseline = 100*(sum(rasterbinned_la(indss,duration),2) - nanmean(sum(rasterbinned_la(indss0,duration),2))) / nanmean(sum(rasterbinned_la(indss0,duration),2));
                                deltaf = (sum(rasterbinned_la(indss,duration),2) - nanmean(sum(rasterbinned_la(indss0,duration),2)));
                                
                                smb_delta(ldelay) = nanmean(deltaf,1);
                                if numel(sum(rasterbinned_la(indss,duration),2))
                                    [smb_p(ldelay),~] =  ranksum(sum(rasterbinned_la(indss0,duration),2),sum(rasterbinned_la(indss,duration),2));
                                end
                                smb_z(ldelay) = (nanmean(sum(rasterbinned_la(indss0,duration),2)) - nanmean(sum(rasterbinned_la(indss,duration),2)))/...
                                    sqrt((nanmean(sum(rasterbinned_la(indss0,duration),2))/length(indss0))+(nanmean(sum(rasterbinned_la(indss,duration),2))/length(indss)));
                                
                                smb(ldelay) = nanmean(differencefromavbaseline,1) ;
                                smb_bs(ldelay) = nanmean(nanmean(sum(rasterbinned_la(indss0,duration),2)),1);
                                smb_ls(ldelay) = nanmean(nanmean(sum(rasterbinned_la(indss,duration),2)),1);
                                %stdmb = std(differencefromavbaseline,0,1);
                                stemb(ldelay) = std(differencefromavbaseline,1)/sqrt(length(indss));
                                stemb_delta(ldelay) = std(deltaf,1)/sqrt(length(indss));
                                
                            end
                            
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%
                            smb_centers(ldelay) = centers(ldelay);
                        end
                    end
                    %%%%% removing Nan entries of smb_centers and corresponding
                    %%%%% for smb, smb_delta, smb_p, smb_bs
                    indsstokill=find(isnan(smb_centers));
                    smb_centers( indsstokill) =[];
                    smb( indsstokill) =[];
                    smb_p( indsstokill) =[];
                    smb_z( indsstokill) =[];
                    smb_bs( indsstokill)  =[];
                    smb_ls( indsstokill)  =[];
                    smb_delta( indsstokill)  =[];
                    nlsAv( indsstokill,:) =[];
                    %laAbs{indsstokill} = [];
                    
                    thiscell.Avfr = Avfr;
                    thiscell.depth = Depth;
                    % setting trialtype as field
                    thiscell.nbs = setfield(thiscell.nbs,trialtype_name,rasterbinnedn(indss0,:));
                    thiscell.nls = setfield(thiscell.nls,trialtype_name,nls);
                    thiscell.nbsAv = setfield(thiscell.nbsAv,trialtype_name,mean(rasterbinnedn(indss0,:),1));
                    thiscell.nlsAv = setfield(thiscell.nlsAv,trialtype_name,nlsAv);
                    % no V1: should be done separately later
                    %thiscell = makecells_makesignoisevorrelation(thiscell);
                    
                    thiscell.smb = setfield(thiscell.smb,trialtype_name,smb);
                    thiscell.smb_p = setfield(thiscell.smb_p,trialtype_name,smb_p);
                    thiscell.smb_z = setfield(thiscell.smb_z,trialtype_name,smb_z);
                    thiscell.smb_bs = setfield(thiscell.smb_bs,trialtype_name,smb_bs);
                    thiscell.smb_ls = setfield(thiscell.smb_ls,trialtype_name,smb_ls);
                    thiscell.smb_delta = setfield(thiscell.smb_delta,trialtype_name,smb_delta);
                    thiscell.smb_centers = setfield(thiscell.smb_centers,trialtype_name,smb_centers);
                    thiscell.averagesilencing = setfield(thiscell.averagesilencing,trialtype_name,nanmean(smb));
                    thiscell.laAbs = setfield(thiscell.laAbs,trialtype_name,laAbs);
                    thiscell.laAls = setfield(thiscell.laAls,trialtype_name,laAls);
                    
                    thisBeh.licks = setfield(thisBeh.licks,trialtype_name,rasterbinnedlicks(indss0,:));
                    thisBeh.licksalltrials = setfield(thisBeh.licksalltrials,trialtype_name,rasterbinnedlicks);
                    thisBeh.ldelayindss = setfield(thisBeh.ldelayindss,trialtype_name,ldelayindss);
                    
                end
                %%%%%%%%%%%%%%%%%%%%%%%%% all thiscell properties are set
                %%%%%%%%%%%%%%%%%%%%%%%%% at this point
                if strcmp(area{:},'LM')
                    LMcells{end+1} = thiscell;
                elseif strcmp(area{:},'V1')
                    V1cells{end+1} = thiscell;
                end
            end
        end
    end
    
    Behcells{end+1} = thisBeh;
end
% saving:
% save('genericcells_ex_40ms_withlicks_groomingonly.mat','LMcells','V1cells','edgestep','Behcells')

% example how to get all V1 Sus for a given LM cell:
% cell2mat(transpose(cellfun(@(x) x(cellfun(@(x) x.simulcode , V1cells)).nbsAv, V1cells,'UniformOutput',0)))