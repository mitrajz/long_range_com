% This function is to assign depth and spike shape to neurons from
% different animals. For depth assignment, make sure you have already done:
%   -LFP_filterdownsample.m : saves time downsampled lfp
%   -LFP_stimulus_induced_LFP_CSD.m : alignes lfp to stimulus in laser and
%    basline trials, calculates csd for basline trials, and estimates depth
%   -LFP_postalignmentchecks_andcorrections: checks depth and csd
%   consistency across animals, possibility or redoing depth assigment.
%
% 
% animallist = {'MPV17','MPV18_2',...
%     'VL53','VL52','VL51','VL66'};%,...
%     %'VL61','VL63','VL55','VL59','VL50'};
% preprocessinglist = {'2019_04_05_21_25_20','2019_04_05_15_04_39',...
% '2017_12_05_18_58_17','2017_12_04_15_08_47','2017_12_03_17_30_15','2018_03_20_11_53_19'};%,...
%    % '2018_03_16_16_53_13','2018_04_01_12_24_38','2018_04_01_15_03_48','2018_04_01_13_24_48','2017_12_01_16_21_49'};
   
%    animallist = {'UMPV2','UMPV4','UMPV7','UMPV8'};
% preprocessinglist = {'2019_09_13_17_56_13','2019_09_13_18_12_17',...
%     '2019_09_13_18_35_50','2019_09_13_19_01_09'};

animallist = {'MPV32_2','MPV33','MPV31','MPV34_2',...
    'MPV36_2','MPV30','MPV35_2'};
preprocessinglist = {'2020_01_21_18_13_31','2020_01_21_17_12_32','2020_01_22_13_47_56','2020_01_22_14_28_52',...
    '2020_01_22_11_38_57','2020_01_22_14_51_15','2020_01_28_20_29_59'};


networkrootpath = '/mnt/javadzam/winstor/swc/mrsic_flogel/public/projects/MiJa_20160601_VisualLongRangeConnectivity/Ephys/';
changepath = 0;
for animal_i=1:length(animallist)
    %% load the preprocessed mat file
    %cd('/mnt/data/Mitra/figs/P2_L/LED/FB');
    cd('/mnt/data/Mitra/figs');
    cd(sprintf('%s/preprocessing/%s',animallist{animal_i},preprocessinglist{animal_i}))
    matfilename = dir('task_*.mat');
    shankfilename = dir('shanki*.mat');
    load(matfilename.name);
    load(shankfilename.name);
    clear valves licks allvalve alllick PStepTimeOn PStepTimeOff LStepTimeOn
    %% extract depth and spikeshapes and append to LM and V1 structs
    areanames = {'V1','LM'};    
    for areai = 1:2
        area = areanames{areai};
        for shank = 1:2    
            % Find the shank depth map
            shankindex = find(strcmp(order,[area,'_',num2str(eval('shank-1'))]));
            aligneddepths = shankdepthmap{shankindex};
            
            % the current electrode:
            TempElectrode = eval(sprintf('%s{shank}',area));       
            % prepare fields for spikeshape and depth, one per each single
            % unit in the shank
            eval(sprintf('%s{shank}.SingleUnitTemplates = cell(1,size(%s{shank}.SingleUnits,1))',area,area));
            eval(sprintf('%s{shank}.SingleUnits_aligned_Depth = nan(1,size(%s{shank}.SingleUnits,1))',area,area));
            % correct path2spikes to new network drive%s{shank}
            old_path = eval('TempElectrode.path2spikes');
            if changepath
                new_path = [networkrootpath,old_path(strfind(old_path,'processed'):end)];
            else
                new_path = old_path;
            end
            % load all spike temples: 
            cd(new_path )
            templatedata= readNPY('templates.npy');
            
            for i = 1:size(TempElectrode.SingleUnits,1)               
                %%% find cluster template id for each cluster :  double
                %%% check the procedure
                unitnum = str2double(TempElectrode.SingleUnits{i}(end-2:end));
                clusterid = ...
                    unique(TempElectrode.tableSpikeData.spike_templates(find(TempElectrode.tableSpikeData.spike_clusters == unitnum)));
                %
                if length(clusterid) > 1
                    numspikes = nan(1,length(clusterid));
                    for clusteri = 1:length(clusterid)
                        potentialclusters = clusterid(clusteri);
                        numspikes(clusteri) = length(find(TempElectrode.tableSpikeData.spike_templates(find(TempElectrode.tableSpikeData.spike_clusters == unitnum)) == potentialclusters));
                    end
                    [~,largerclusterid]=max(numspikes);
                    clusterid = clusterid(largerclusterid);
                end
                % this is time*space; profile of spike over channels
                st_template = squeeze(templatedata(clusterid+1,:,:));                
                % find the largest magnitute channel
                [~,ch]=max(max(abs(st_template)));
                
                % assign spikeshape to the field
                eval(sprintf('%s{shank}.SingleUnitTemplates{i} = st_template(:,ch);',area))
                % assign depth
                eval(sprintf('%s{shank}.SingleUnits_aligned_Depth(i) = aligneddepths(TempElectrode.SingleUnitsDepths(i));',area))
            end
            
            figure;plot(cell2mat(eval(sprintf('%s{shank}.SingleUnitTemplates',area))))
            title(sprintf('%s shank%d ',area,shank))
        end
    end
    %% replace V1 and LM structs in preprocessed mat file   
    cd('/mnt/data/Mitra/figs');
    cd(sprintf('%s/preprocessing/%s',animallist{animal_i},preprocessinglist{animal_i}))
    save(matfilename.name,'LM','V1','-append','-v7.3')
end
