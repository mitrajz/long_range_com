% Before running this script, kilosort retin recordings (no manual sorting
% needed), and then run spikeshape_comparison.m 
% The script goes over all
% animals, takes task file details, and retin conversion, and makes and
% saves retin_area{shank} (saves as retinmaps.m in the task
% preprocessing folder) , and retinotopy plots for each animals (saved in
% retin_figures)
% No time stamp, if files already exist, they get overwritten
% retin maps are done with 2D gaussian fits, default parameters are on and
% off fields averaged.
%%
clear all
%%
addpath(genpath('/mnt/data/Mitra/cache/codes/matlabkwikreader/'))
addpath(genpath('/mnt/data/Mitra/cache/codes/'))
%%
animallist = {'MPV17','MPV18_2',...
    'VL53','VL52','VL51','VL66',...
    'VL61','VL63','VL55','VL59','VL50',...
    'MPV32_2','MPV33','MPV31','MPV34_2',...
    'MPV36_2','MPV30','MPV35_2'};
preprocessinglist = {'2019_04_05_21_25_20','2019_04_05_15_04_39',...
    '2017_12_05_18_58_17','2017_12_04_15_08_47','2017_12_03_17_30_15','2018_03_20_11_53_19',...
    '2018_03_16_16_53_13','2018_04_01_12_24_38','2018_04_01_15_03_48','2018_04_01_13_24_48','2017_12_01_16_21_49',...
    '2020_01_21_18_13_31','2020_01_21_17_12_32','2020_01_22_13_47_56','2020_01_22_14_28_52',...
    '2020_01_22_11_38_57','2020_01_22_14_51_15','2020_01_28_20_29_59'};

processed_root = '/mnt/data/Mitra/processed/';
raw_root = '/mnt/data/Mitra/raw/';

gfit.on_off = nan; % on_off: 0:off, 1:on, nan: combine
gfit.sres = 0.01;% spatial resolution for evaluating gaussian fits.
gfit.ResponseOnset_ms = 50; % Response to the stimulus is calculated from this point to end of Window which is 200ms default
Normalize = 0; % normalize each units map to its max or not, just for plotting. Doesn't effect saved maps.

for animal_i=1:length(animallist)
    
    animallist{animal_i}
    
    cd('/mnt/data/Mitra/figs');
    cd(sprintf('%s/preprocessing/%s',animallist{animal_i},preprocessinglist{animal_i}))
    
    % if conversion map already exists, ask if you want to re-do
    if ~exist('tasktoretin_unit_conv.mat','file')
        disp('no conversion map');break;
    end
    
    
    matfilename = dir('task_*.mat');
    load(matfilename.name)
    retinconvfilename = dir('tasktoretin_unit_conv.mat');
    load(retinconvfilename.name)
    
    clear tStart tEnd
    WindowSec = .2;
    % Retinpath is only for 1 shank, cut up and do for all
    % construct paths:    
    splt = strsplit(Retinpath,'/');
    splt2 = strsplit(splt{end},'_kwd');
    R_RootPath = splt2{1};
    
    % be careful that V1,LM already exist (task loaded), so use prefix retin
    retin_kwd = [raw_root,animallist{animal_i},'/',R_RootPath,'/experiment1_101.raw.kwd'];
    retin_V1{1}.path2spikes = [processed_root,animallist{animal_i},'/',R_RootPath,'_kwd_101_rec0_shkV1_shk0/'];
    retin_V1{2}.path2spikes =  [processed_root,animallist{animal_i},'/',R_RootPath,'_kwd_101_rec0_shkV1_shk1/'];
    retin_LM{1}.path2spikes = [processed_root,animallist{animal_i},'/',R_RootPath,'_kwd_101_rec0_shkLM_shk0/'];
    retin_LM{2}.path2spikes =  [processed_root,animallist{animal_i},'/',R_RootPath,'_kwd_101_rec0_shkLM_shk1/'];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % load raw data, read pd and make PDAllOn
    
    
    PdTrig_nofilter = h5read(retin_kwd,kwdrecording,[Pdch 1],[1 inf]);
    figure;plot(PdTrig_nofilter(1:1*10^6));
    userin = input('tStart:','s');
    tStart = str2num(userin)
    
    figure;plot(PdTrig_nofilter(end-10^5:end));
    userin = input('How much to cut:','s');
    tEnd = length(PdTrig_nofilter) - str2num(userin)
    
    
    filt_pd=low_pass_filter(PdTrig_nofilter,samplingf,pdfiltfreq);
    [StepTimeOn,StepTimeOff,StepTimeStampOn,StepTimeStampOff] = FindEdges(filt_pd,DiffTh,tStart,tEnd);
    
    sprintf('%d pd onset detected',length(StepTimeStampOn))
    sprintf('%d pd offset detected',length(StepTimeStampOff))
    
    Window=WindowSec*samplingf;
    PDAllOn= repmat(StepTimeStampOn',[1,2*Window+1])+repmat(double(-Window:Window),[length(StepTimeStampOn'),1]); % if you don't double, it will be single and the subtraction would be fucked
    
    figure;plot((-Window:Window)/samplingf,filt_pd(PDAllOn(1:40,:))');
    
    % 6 or 8 repetitions of 40 ordered stimuli, one round white, other black
    if length(StepTimeStampOn) == 240
        Order =  repmat(1:40,[1,6]);
    elseif length(StepTimeStampOn) == 320
        Order =  repmat(1:40,[1,8]);
    else
        disp('something wrong with pd detection');break;
    end
    %%%%%%%%%%%%%
    % Read spikes and make retin maps using gaussian fit
    %%%% for an example shank
    targetareanames = {'V1','LM'};
    
    for ti = 1:length(targetareanames)
        targetarea = eval(['retin_',targetareanames{ti}]);
        unit_maps_sh = cell(1,2);
        
        for shank = 1:2
            AllSpikes=[readNPY([targetarea{shank}.path2spikes,'spike_times.npy'])];
            clusters=readNPY([targetarea{shank}.path2spikes, 'spike_clusters.npy']); 
            
            % main function performing gaussian fit and assigning
            % retinotopy maps:       
            unit_maps_sh{shank} = retin_process_unit(AllSpikes,clusters,Order,Window,PDAllOn,filt_pd,gfit.on_off,gfit.sres,gfit.ResponseOnset_ms);
            
            % test for single units
            if 0
                cluster = 3;
                figure;subplot(2,1,1);imagesc(unit_maps.RespMag{cluster})
                colormap('gray')
                subplot(2,1,2);im=imagesc([0.5,5.5],[0.5,4.5],unit_maps.predictedmap{cluster});
                im.AlphaData = .5;
                axis tight
                
                hold on;plot(unit_maps.xgcenter{cluster},unit_maps.ygcenter{cluster},'r+')
            end
            
            %%%%%%%%%%% making indices for units
            try % if cluster_KSLabel.tsv exists
                KSL=tdfread([targetarea{shank}.path2spikes,'cluster_KSLabel.tsv']);
                %%% good KS labels:  check first 2 letters
                goodclsters_ind = KSL.cluster_id(find((arrayfun(@(x) strcmp(x,'g'),KSL.KSLabel(:,1)) & arrayfun(@(x) strcmp(x,'o'),KSL.KSLabel(:,2)))));
            catch % if cluster_KSLabel.tsv doesn't exist
                goodclsters_ind = [];
            end
            
            %%% inices of matching to task
            fieldname = [targetareanames{ti},'_',num2str(shank-1)];
            matchedunit_ind = getfield(task_retin_SU_conversion,fieldname);
            % conversion indices are indices of templates (not name of unit), use the line bellow to
            % convert them to relevent format (element numbet of unit_maps_sh{shank} fields)
            % matchedunit_ind should have as many elemnts as the task
            % shank, with each element directly indexing unit_maps_sh{shank}.clusterid
            % LM{1}.sigleunit{2} -> retin info, e.g.xg center in:
            % unit_maps_sh{shank}.xgcenter{matchedunit_ind(2)},...
            matchedunit_ind = convert_tableunit_to_templateunit(matchedunit_ind,targetarea{shank}.path2spikes);     
            % For here, the exact matchinf doesn't matter. So:
            matchedunit_ind_nonan = matchedunit_ind(find(~isnan(matchedunit_ind)));
                        
            %%% assign indices as fields
            unit_maps_sh{shank}.matchedunit_ind = matchedunit_ind;
            unit_maps_sh{shank}.goodclsters_ind = goodclsters_ind';
            
            %%%% making flags:
            %%%% Choosing units: a binary field to include the unit or not
            unit_maps_sh{shank}.have_retin_flag = cellfun(@(x) numel(x),unit_maps_sh{shank}.predictedmap)>1;            
            unit_maps_sh{shank}.goodKS_flag = cellfun(@(x) numel(intersect(x,unit_maps_sh{shank}.goodclsters_ind)),unit_maps_sh{shank}.clusterid);
            unit_maps_sh{shank}.matchedunit_flag = cellfun(@(x) numel(intersect(x,matchedunit_ind_nonan)),unit_maps_sh{shank}.clusterid);

            % save in the corresponding retin_area{shank} field
            eval(sprintf('retin_%s{shank}.retinmap = unit_maps_sh{shank}',targetareanames{ti}))            
        end                
    end
    
    
    [retin_V1,retin_LM] = calculatedistances(retin_V1,retin_LM);
    
    
    %%% merge 2 shanks of each area and make figures
    f_m=figure;
    mergeshanksandplot({retin_V1{1}.retinmap,retin_V1{2}.retinmap},'taskmatched',subplot(2,2,1),subplot(2,2,3),'V1',Normalize)
    mergeshanksandplot({retin_LM{1}.retinmap,retin_LM{2}.retinmap},'taskmatched',subplot(2,2,2),subplot(2,2,4),'LM',Normalize)
    f_m.Name = [animallist{animal_i},'-','taskmatched'];
    
    f_g=figure;
    mergeshanksandplot({retin_V1{1}.retinmap,retin_V1{2}.retinmap},'KSgood',subplot(2,2,1),subplot(2,2,3),'V1',Normalize)
    mergeshanksandplot({retin_LM{1}.retinmap,retin_LM{2}.retinmap},'KSgood',subplot(2,2,2),subplot(2,2,4),'LM',Normalize)
    f_g.Name = [animallist{animal_i},'-','KSgood'];

    f_all=figure;
    mergeshanksandplot({retin_V1{1}.retinmap,retin_V1{2}.retinmap},'all',subplot(2,2,1),subplot(2,2,3),'V1',Normalize)
    mergeshanksandplot({retin_LM{1}.retinmap,retin_LM{2}.retinmap},'all',subplot(2,2,2),subplot(2,2,4),'LM',Normalize)
    f_all.Name = [animallist{animal_i},'-','all'];

    % save retin_area and figure handles, maybe make a retin folder
    %%% save
    cd('/mnt/data/Mitra/figs');
    cd(sprintf('%s/preprocessing/%s',animallist{animal_i},preprocessinglist{animal_i}))
    save('retinmaps.mat','retin_V1','retin_LM')
    % save retin figures
    if ~exist('retin_figures','dir')
        mkdir('retin_figures')
    end
    cd('retin_figures')
    
    hgsave(f_m,['retin_',f_m.Name,'.fig'], '-v7.3');
    saveas(f_m,['retin_',f_m.Name,'.png']);
    hgsave(f_g,['retin_',f_g.Name,'.fig'], '-v7.3');
    saveas(f_g,['retin_',f_g.Name,'.png']);
    hgsave(f_all,['retin_',f_all.Name,'.fig'], '-v7.3');
    saveas(f_all,['retin_',f_all.Name,'.png']);
       
end

