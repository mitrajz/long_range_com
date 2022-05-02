% IMPORTANT: OUTPUT: conversion file is the INDEX of retin_template which
% consists of the names of all units that have spikes. It is different from
% spike template and KS label for example. When loading the conversion use 
% convert_tableunit_to_templateunit to convertbetween different indices
%
% Matching single units of task recording to the independently sorted retin
% units.
% Before running: manual sorting of retin

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

spikes_rootfile = '/mnt/data/Mitra/processed/';

for animal_i=1:length(animallist)
    animallist{animal_i}
    
    
    
    proceed = 1;
    
    cd('/mnt/data/Mitra/figs');
    cd(sprintf('%s/preprocessing/%s',animallist{animal_i},preprocessinglist{animal_i}))
    
    % if conversion map already exists, ask if you want to re-do
    if exist('tasktoretin_unit_conv.mat','file')
        proceed = input('redo?[0/1]');
    end
    
    if proceed
        
        matfilename = dir('task_*.mat');
        load(matfilename.name)
        
        % finding task and retin file names
        Lfilename = [expname{1}(2:end),'_kwd_101_rec0_'];
        expname{1}(2:end)
        userretinfolder = uigetdir(['/mnt/data/Mitra/raw/',animallist{animal_i}]);
        % last bit of the name
        seps = strfind(userretinfolder,'/');
        rawretinname = userretinfolder(seps(end)+1:end);
        Retinfilename = [rawretinname,'_kwd_101_rec0_'];
        
        task_retin_SU_conversion = struct;
        
        areanames = {'V1','LM'};
        for areanum = 1:2
            for shanknum = 1:2
                compound_name = [areanames{areanum},'_',num2str(shanknum-1)]
                Lpath  = [spikes_rootfile,animallist{animal_i},'/',Lfilename,'shk',areanames{areanum},'_shk',num2str(shanknum-1)];
                Retinpath  = [spikes_rootfile,animallist{animal_i},'/',Retinfilename,'shk',areanames{areanum},'_shk',num2str(shanknum-1)];
                              
                % template for good spikes in task and retin recordings,
                % unit_type = either 'good' or 'unsorted'
                L_st_template = getspiketemplate_for_goodunits(Lpath,'good'); % unit_type, either 'good or 'all'
                Retin_st_template = getspiketemplate_for_goodunits(Retinpath,'all');
                
                
                
                % array of length = number of single units in task recording. each unit
                % gives the correspnding index in retin recording.
                this_task_retin_SU_conversion = nan(1,length(L_st_template));
                try
                    for n= 1:length(L_st_template)
                        % difference between unit n from task to all units in retin
                        distance = nan(1,length(Retin_st_template));
                        % different distance measures: watch out with
                        % precision! Variables are single
                        for r = 1:length(Retin_st_template)
                           
                            distance(r) = nansum(nansum((L_st_template{n} - Retin_st_template{r}).^2));
                          
                        end
                        % MSE of match to all retin single cells
                        [~,mostsimilarind] = min(distance);
                        [~,similarityind]=sort(distance);
                        % manual inspection. If not approved, index remains nan
                        % 3 trials to get to best match
                        for u = 1:3
                            testf = figure;
                            plot(L_st_template{n},'k');hold on;plot(Retin_st_template{similarityind(u)},'r');
                            manualmatch = input('accept?[0/1]');
                            if manualmatch
                                this_task_retin_SU_conversion(n) = similarityind(u);%mostsimilarind;
                                close(testf);
                                break;
                            else
                                close(testf);
                            end
                            
                        end
                        
                    end
                catch
                    disp('didnt make conversion file for this shank') 
                end
                task_retin_SU_conversion = setfield(task_retin_SU_conversion,compound_name,this_task_retin_SU_conversion);
            end
        end
        
        
        
        % go to preprocessing folder
        cd('/mnt/data/Mitra/figs');
        cd(sprintf('%s/preprocessing/%s',animallist{animal_i},preprocessinglist{animal_i}))
        % save file
        save('tasktoretin_unit_conv.mat','task_retin_SU_conversion','Retinpath');
    end
end

function st_template = getspiketemplate_for_goodunits(Lpath,unit_type)

[tableSpikeData, tableClusterData, strParams] = loadKiloSortedSpikes(Lpath);



ClusterTypes = fieldnames(tableClusterData);
AllUnits=[];
SingleUnits=[];
tableClusterData.all=[];
AllUnitSpikeTimes=cell(1,length(AllUnits));
SingleUnitSpikeTimes=cell(1,length(AllUnits));
if strcmp(unit_type,'good')
    for i=1:length(ClusterTypes)
        if strcmp(ClusterTypes{i},'good')
            SingleUnits = [SingleUnits ; fieldnames(getfield(tableClusterData,ClusterTypes{i}))];
            ClusterUnitSpikes = cellfun(@(x) getfield(getfield(tableClusterData,ClusterTypes{i}),x) , fieldnames(getfield(tableClusterData,ClusterTypes{i})),'UniformOutput', false);
            SingleUnitSpikeTimes = [SingleUnitSpikeTimes, ClusterUnitSpikes'];
        end
    end
else
    for i=1:length(ClusterTypes)
        SingleUnits = [SingleUnits ; fieldnames(getfield(tableClusterData,ClusterTypes{i}))];
        ClusterUnitSpikes = cellfun(@(x) getfield(getfield(tableClusterData,ClusterTypes{i}),x) , fieldnames(getfield(tableClusterData,ClusterTypes{i})),'UniformOutput', false);
        SingleUnitSpikeTimes = [SingleUnitSpikeTimes, ClusterUnitSpikes'];
    end
end
%%%

templatedata= readNPY(fullfile(Lpath,'templates.npy'));
st_template = cell(1,size(SingleUnits,1));

for i = 1:size(SingleUnits,1)
    %%% find cluster template id for each cluster :  double
    %%% check the procedure
    unitnum = str2double(SingleUnits{i}(end-2:end));
    clusterid = ...
        unique(tableSpikeData.spike_templates(find(tableSpikeData.spike_clusters == unitnum)));

    if length(clusterid) > 1
        numspikes = nan(1,length(clusterid));
        for clusteri = 1:length(clusterid)
            potentialclusters = clusterid(clusteri);
            numspikes(clusteri) = length(find(tableSpikeData.spike_templates(find(tableSpikeData.spike_clusters == unitnum)) == potentialclusters));
        end
        [~,largerclusterid]=max(numspikes);
        clusterid = clusterid(largerclusterid);
    end
    % this is time*space; profile of spike over channels
    st_template{i} = squeeze(templatedata(clusterid+1,:,:));
   
end
end