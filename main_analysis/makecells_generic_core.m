function [LMcells,V1cells,Behcells] = makecells_generic_core(params,animallist,preprocessinglist,exptype,trialtypelist,targettrialtypelist)
% make generic LM and V1 cells (both LM and V1), static mice (running speed = 0)

LMcells = cell(0,0);
V1cells = cell(0,0);
Behcells = cell(0,0);



% quick performance check: 
% first lick
%[size(intersect(misstrialind,nogroomingind))/size(intersect(gotrialind,nogroomingind)) size(intersect(falsetrialind,nogroomingind))/size(intersect(nogotrialind,nogroomingind))  ]

trialtype_namelist = {'go','nogo'}; % determines the name of struct fields
for animal_i=1:length(animallist) 
    cd('/mnt/data/Mitra/figs');
    cd(sprintf('%s/preprocessing/%s',animallist{animal_i},preprocessinglist{animal_i}))
    matfilename = dir('task_*.mat');
    perfmatfilename = dir('rev_perf*.mat');
    retinfilename = dir('retinmaps.mat');
    if ~isempty(perfmatfilename)
        load(perfmatfilename.name)
    end
    if ~isempty(retinfilename)
        load(retinfilename.name)
    end
    load(matfilename.name)
    
    if params.Window ~= Window
        error('window size mismatch')
    end
    %%%%%%%%%%%%%%%%%%%%%%%% IMPRTANT: define params after loading mat
    %%%%%%%%%%%%%%%%%%%%%%%% files (to overwrite)
    % assing params from params input
    %clear stabletrialind
    edgestep_pl = params.edgestep_pl;
    alignstyle_pl = params.alignstyle_pl;
    edgestep_an = params.edgestep_an;
    min_n_trials_per_delay = params.min_n_trials_per_delay;
    %
    excludegrooming = params.excludegrooming;
    onlycorrect = params.onlycorrect;
    onlystableperiod = params.onlystableperiod;
    %
    movingwin_withlags.lagsize = params.movingwin_withlags.lagsize;
    movingwin_withlags.nlags = params.movingwin_withlags.nlags;
    movingwin_withlags.movingwinsize = params.movingwin_withlags.movingwinsize;
    %
    if isfield(params,'onsetlag_an') && ~isnan(params.onsetlag_an)
        onsetlag_an = params.onsetlag_an;
    else
        onsetlag_an = 0;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%
    for area = {'V1','LM'}
        for shank=1:2
            for neuronnum =1:size(eval(sprintf('%s{shank}.SingleUnitSpikeTimes',area{:})),2)
                %%%%%%%%%%%%%%%%%%%%%%%%%% prepare empty struct fields
                %%%%%% class 1 : structural properties of the cell,
                %%%%%% independent of most things
                clear thiscell
                % trialtype independent
                thiscell.animal = animalname;
                thiscell.simulcode = animal_i;
                thiscell.exptype = exptype{animal_i};
                thiscell.shank = shank;
                thiscell.neuronnum = neuronnum;
                thiscell.expname = expname;
                thiscell.depth = eval(sprintf('%s{shank}.SingleUnits_aligned_Depth(neuronnum)',area{:}));
                thiscell.spiketemplate = eval(sprintf('%s{shank}.SingleUnitTemplates{neuronnum}',area{:}));

                thiscell.retin = struct;
                %%%%% class 1.5: mostly independent
                thiscell.stimspikes =  struct; % spikes during stim (500ms)
                thiscell.prestimspikes =  struct; % spikes pre stim (500ms)
                %%%%%% class 2: trialtype dependent - plotting fields -
                %%%%%% not to be used for analysis
                thiscell.nbs = struct;
                thiscell.nls = struct;
                thiscell.nlsAv = struct;
                thiscell.nbsAv = struct;
                %%%%%% class 3: trialtype dependent: analysis fields - 
                %%%%%% characterizing silencing
                thiscell.smb = struct;
                thiscell.smb_ctrl = struct;
                thiscell.smb_p = struct;
                thiscell.smb_z = struct;
                thiscell.smb_bs = struct;
                thiscell.smb_ls = struct;
                thiscell.smb_delta = struct;
                thiscell.smb_centers = struct;
                thiscell.laAbs = struct;
                thiscell.laAls = struct;
                thiscell.laAbs_ctrl = struct;
                thiscell.laAls_ctrl = struct;
                thiscell.laAcontls = struct;
                %%%%%% class 4: lag analysis fields: used for onset-offset 
                thiscell.effectsize_lags = struct;
                thiscell.effectsize_lags_pre = struct;
                thiscell.pval_lags = struct;
                %%%%%% licks use plotting bin - should not be used for
                %%%%%% analysis in the current form
                thisBeh.licks = struct;
                thisBeh.licksalltrials = struct;
                thisBeh.ldelayindss = struct;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%% set up retin fields
                targetretin = eval(['retin_',area{:}]);
                % find the matched index
                clusterid_index = targetretin{shank}.retinmap.matchedunit_ind(neuronnum);
                retin_index = find(cellfun(@(x) x == clusterid_index,targetretin{shank}.retinmap.clusterid));
                % if match exists, and so does a gaussian fit for the cell add center, local, and longrange centers
                if numel(retin_index) && (numel(targetretin{shank}.retinmap.xgcenter{retin_index})>0)
                    
                    thiscell.retin.center = [targetretin{shank}.retinmap.xgcenter{retin_index},...
                        targetretin{shank}.retinmap.ygcenter{retin_index}];
                    thiscell.retin.distancelocal = targetretin{shank}.retinmap.distancelocal.distance{retin_index};
                    thiscell.retin.distancelong = targetretin{shank}.retinmap.distancelong.distance{retin_index};
                    
                % otherwise fill the fields with nan
                else
                    thiscell.retin.center = [nan nan];
                    thiscell.retin.distancelocal = nan;
                    thiscell.retin.distancelong = nan;
                end
                %%%%%%%%
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
                    Avfr = samplingf*length(units)/double(units(end) - units(1));
                    stimdur = 500*30;
                    
                    if ~(alignstyle_pl) % old style, align to -Window
                        rasterbinnedn=nan(size(PAllOn,1),(length(1:edgestep_pl:size(PAllOn,2))-1));
                    elseif alignstyle_pl % new style, align to stim onset
                        % an example to set the size of rasterbinned
                        edges_size = [fliplr(PAllOn(1,ceil(size(PAllOn,2)/2))-edgestep_pl:-edgestep_pl:PAllOn(1,1)),PAllOn(1,ceil(size(PAllOn,2)/2)):edgestep_pl:PAllOn(1,end)];
                        rasterbinnedn=nan(size(PAllOn,1),length(edges_size)-1);
                    end
                    
                    rasterbinned_allstim = nan(size(PAllOn,1),2);
                    rasterbinnedlicks=nan(size(PAllOn,1),(length(1:edgestep_pl:size(PAllOn,2))-1));
                    for i=1:1:size(PAllOn,1)
                        middlep = ceil(size(PAllOn,2)/2);
                        if ~(alignstyle_pl)
                            edges = PAllOn(i,1):edgestep_pl:PAllOn(i,end);
                        elseif alignstyle_pl
                            %%%% IMPORTANT: Window is assumed symetric here
                            % ceil or floor makes 1/30000 sec of
                            % difference, it is not important
                            edges = [fliplr(PAllOn(i,middlep-edgestep_pl):-edgestep_pl:PAllOn(i,1)),PAllOn(i,middlep):edgestep_pl:PAllOn(i,end)];
                        end
                        [rasterbinnedn(i,:),~] = histcounts(units,edges);
                        
                        [rasterbinned_allstim(i,1),~] = histcounts(units,PAllOn(i,middlep)-stimdur:stimdur:PAllOn(i,middlep));
                        [rasterbinned_allstim(i,2),~] = histcounts(units,PAllOn(i,middlep):stimdur:PAllOn(i,middlep)+stimdur);
                        
                        edgeslick = 0:edgestep_pl:size(licks,2)-1;
                        [rasterbinnedlicks(i,:),~] = histcounts(find(licks(i,:)),edgeslick);
                    end

       
                    %%%%%% silencing effect
                    clear smb stemb
                    smb = nan(1,max(unique(LaserDelayBinned)));
                    smb_ctrl = nan(1,max(unique(LaserDelayBinned)));
                    smb_p = nan(1,max(unique(LaserDelayBinned)));
                    smb_z = nan(1,max(unique(LaserDelayBinned)));
                    smb_bs = nan(1,max(unique(LaserDelayBinned)));
                    smb_ls = nan(1,max(unique(LaserDelayBinned)));
                    stemb = nan(1,max(unique(LaserDelayBinned)));
                    smb_centers = nan(1,max(unique(LaserDelayBinned)));
                    laAbs = cell(0,0);
                    laAls = cell(0,0);
                    laAcontls = cell(0,0);
                    laAbs_ctrl = cell(0,0);
                    laAls_ctrl = cell(0,0);
                    nls = cell(0,0);
                    effectsize_lags = cell(0,0);
                    effectsize_lags_pre = cell(0,0);
                    pval_lags = cell(0,0);
                    
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
                            nls{end+1} = rasterbinnedn(indss,:);
                            %%%%%%%%%%%%%%%%%% 
                            % laser aligned activity:
                            % activity in 1 bin after laser onset, in both
                            % baseline and laser trials (except for cont
                            % which contains mutiple bins)
                            rasterbinned_la=nan(size(PAllOn,1),1);
                            rasterbinned_la_ctrl=nan(size(PAllOn,1),1);
                            rasterbinned_la_cont=nan(size(PAllOn,1),ceil(stimdur/edgestep_an));
                            indsstofill = [indss0;indss];
                            point=floor(floor(size(PAllOn,2)/2)+centers(ldelay)*30)+onsetlag_an;
                            for i=1:1:size(indsstofill,1)
                                edges =PAllOn(indsstofill(i),point):edgestep_an:PAllOn(indsstofill(i),point+1*edgestep_an);
                                edges_ctrl =PAllOn(indsstofill(i),point-1*edgestep_an):edgestep_an:PAllOn(indsstofill(i),point);
                                edges_con = PAllOn(indsstofill(i),point):edgestep_an:...
                                    PAllOn(indsstofill(i),floor(floor(size(PAllOn,2)/2)+stimdur));
                                
                                [rasterbinned_la(indsstofill(i),:),~] = histcounts(units,edges);
                                [rasterbinned_la_ctrl(indsstofill(i),:),~] = histcounts(units,edges_ctrl);
                                if length(edges_con)>1
                                [rasterbinned_la_cont(indsstofill(i),1:length(histcounts(units,edges_con))),~] = histcounts(units,edges_con);
                                end
                            end
                            laAbs{end+1} = rasterbinned_la(indss0);
                            laAls{end+1} = rasterbinned_la(indss);
                            laAcontls{end+1} = rasterbinned_la_cont(indss,:);
                            laAbs_ctrl{end+1} = rasterbinned_la_ctrl(indss0);
                            laAls_ctrl{end+1} = rasterbinned_la_ctrl(indss);
                            %%%%%%%%%%%%%%%%% 
                            % detecting the onset of effect after laser: a
                            % window is moved with different lags from
                            % laser onset, to see which lag maximizes the
                            % effect size of laser and minimizes the
                            % pvalue. parameters are set in movingwin_withlags
                            
                            rasterbinned_lags=nan(size(PAllOn,1),movingwin_withlags.nlags); 
                            rasterbinned_lags_pre=nan(size(PAllOn,1),movingwin_withlags.nlags); 
                            indsstofill = [indss0;indss];
                            point=floor(floor(size(PAllOn,2)/2)+centers(ldelay)*30);
                            for i=1:1:size(indsstofill,1)
                                for startlag = 1:movingwin_withlags.nlags
                                    edges =PAllOn(indsstofill(i),point+(startlag-1)*movingwin_withlags.lagsize):movingwin_withlags.movingwinsize:PAllOn(indsstofill(i),point+(startlag-1)*movingwin_withlags.lagsize+movingwin_withlags.movingwinsize);
                                    [rasterbinned_lags(indsstofill(i),startlag),~] = histcounts(units,edges);
                                end
                                for startlag = 1:movingwin_withlags.nlags
                                    edges =PAllOn(indsstofill(i),point-(startlag-1)*movingwin_withlags.lagsize-movingwin_withlags.movingwinsize):movingwin_withlags.movingwinsize:PAllOn(indsstofill(i),point-(startlag-1)*movingwin_withlags.lagsize);
                                    [rasterbinned_lags_pre(indsstofill(i),startlag),~] = histcounts(units,edges);
                                end
                            end
                            effectsize_lags{end+1} = arrayfun(@(x) (nanmean(rasterbinned_lags(indss,x))-nanmean(rasterbinned_lags(indss0,x))),1:movingwin_withlags.nlags);
                            effectsize_lags_pre{end+1} = arrayfun(@(x) (nanmean(rasterbinned_lags_pre(indss,x))-nanmean(rasterbinned_lags_pre(indss0,x))),1:movingwin_withlags.nlags);
                            pval_lags{end+1} = arrayfun(@(x) ranksum(rasterbinned_lags(indss0,x),rasterbinned_lags(indss,x)),1:movingwin_withlags.nlags);
                            
                            %%%%%%%%%%%%%%%%% 
                            % quantifying the effect of laser with laser
                            % aligned activity in laser and baseline trials
                            duration = 1;
                                                        
                            differencefromavbaseline = 100*(rasterbinned_la(indss,:) - nanmean(rasterbinned_la(indss0,:),1)) / nanmean(rasterbinned_la(indss0,:),1);
                            differencefromavbaseline_ctrl = 100*(rasterbinned_la_ctrl(indss,:) - nanmean(rasterbinned_la_ctrl(indss0,:),1)) / nanmean(rasterbinned_la_ctrl(indss0,:),1);
                            deltaf = (rasterbinned_la(indss,:) - nanmean(rasterbinned_la(indss0,:),1));
                            
                            smb_delta(ldelay) = nanmean(deltaf,1);
                            if numel(rasterbinned_la(indss,:)) && numel(rasterbinned_la(indss0,:))
                                [smb_p(ldelay),~] =  ranksum(rasterbinned_la(indss0,:),rasterbinned_la(indss,:));
                            end
                            % double check 
                            smb_z(ldelay) = (nanmean(rasterbinned_la(indss0,:)) - nanmean(rasterbinned_la(indss,:)))/...
                                sqrt((nanmean(rasterbinned_la(indss0,:))/length(indss0))+(nanmean(rasterbinned_la(indss,:))/length(indss)));
                            
                            smb(ldelay) = nanmean(differencefromavbaseline,1) ;
                            smb_ctrl(ldelay) = nanmean(differencefromavbaseline_ctrl,1) ;
                            smb_bs(ldelay) = nanmean(rasterbinned_la(indss0,:),1);
                            smb_ls(ldelay) = nanmean(rasterbinned_la(indss,:),1);
                            %stdmb = std(differencefromavbaseline,0,1);
                            stemb(ldelay) = std(differencefromavbaseline,1)/sqrt(length(indss));
                            stemb_delta(ldelay) = std(deltaf,1)/sqrt(length(indss));
                            %%%%%%%%%%%%%%%%%
                            smb_centers(ldelay) = centers(ldelay);
                        end
                    end
                    %%%%% removing Nan entries of smb_centers and corresponding
                    %%%%% for smb, smb_delta, smb_p, smb_bs
                    indsstokill=find(isnan(smb_centers));
                    smb_centers( indsstokill) =[];
                    smb( indsstokill) =[];
                    smb_ctrl( indsstokill) =[];
                    smb_p( indsstokill) =[];
                    smb_z( indsstokill) =[];
                    smb_bs( indsstokill)  =[];
                    smb_ls( indsstokill)  =[];
                    smb_delta( indsstokill)  =[];
                    nlsAv( indsstokill,:) =[];
                    %laAbs{indsstokill} = [];
                    
                    thiscell.Avfr = Avfr;                 
                    thiscell.prestimspikes =  setfield(thiscell.prestimspikes,trialtype_name,rasterbinned_allstim(indss0,1));
                    thiscell.stimspikes =  setfield(thiscell.stimspikes,trialtype_name,rasterbinned_allstim(indss0,2));
                    % setting trialtype as field
                    thiscell.nbs = setfield(thiscell.nbs,trialtype_name,rasterbinnedn(indss0,:));
                    thiscell.nls = setfield(thiscell.nls,trialtype_name,nls);
                    thiscell.nbsAv = setfield(thiscell.nbsAv,trialtype_name,mean(rasterbinnedn(indss0,:),1));
                    thiscell.nlsAv = setfield(thiscell.nlsAv,trialtype_name,nlsAv);
                    % no V1: should be done separately later
                    %thiscell = makecells_makesignoisevorrelation(thiscell);
                    
                    thiscell.smb = setfield(thiscell.smb,trialtype_name,smb);
                    thiscell.smb_ctrl = setfield(thiscell.smb_ctrl,trialtype_name,smb_ctrl);
                    thiscell.smb_p = setfield(thiscell.smb_p,trialtype_name,smb_p);
                    thiscell.smb_z = setfield(thiscell.smb_z,trialtype_name,smb_z);
                    thiscell.smb_bs = setfield(thiscell.smb_bs,trialtype_name,smb_bs);
                    thiscell.smb_ls = setfield(thiscell.smb_ls,trialtype_name,smb_ls);
                    thiscell.smb_delta = setfield(thiscell.smb_delta,trialtype_name,smb_delta);
                    thiscell.smb_centers = setfield(thiscell.smb_centers,trialtype_name,smb_centers);
                    thiscell.laAbs = setfield(thiscell.laAbs,trialtype_name,laAbs);
                    thiscell.laAls = setfield(thiscell.laAls,trialtype_name,laAls);
                    thiscell.laAbs_ctrl = setfield(thiscell.laAbs_ctrl,trialtype_name,laAbs_ctrl);
                    thiscell.laAls_ctrl = setfield(thiscell.laAls_ctrl,trialtype_name,laAls_ctrl);
                    thiscell.laAcontls = setfield(thiscell.laAcontls,trialtype_name,laAcontls);
                    
                    thiscell.effectsize_lags = setfield(thiscell.effectsize_lags,trialtype_name,effectsize_lags);
                    thiscell.effectsize_lags_pre = setfield(thiscell.effectsize_lags_pre,trialtype_name,effectsize_lags_pre);
                    thiscell.pval_lags = setfield(thiscell.pval_lags,trialtype_name,pval_lags);
                    
                    
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