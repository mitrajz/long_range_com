function [V1,LM] = Extractspikeinfo2(V1,LM)
%%% LM
for shank = 1:length(LM)
    [LM{shank}.tableSpikeData, LM{shank}.tableClusterData, LM{shank}.strParams] = loadKiloSortedSpikes(LM{shank}.path2spikes);
    ClusterTypes = fieldnames(LM{shank}.tableClusterData);
    LM{shank}.AllUnits=[];
    LM{shank}.SingleUnits=[];
    LM{shank}.tableClusterData.all=[];
    LM{shank}.AllUnitSpikeTimes=cell(1,length(LM{shank}.AllUnits));
    LM{shank}.SingleUnitSpikeTimes=cell(1,length(LM{shank}.AllUnits));
    for i=1:length(ClusterTypes)
        if ~strcmp(ClusterTypes{i},'noise')
            LM{shank}.AllUnits = [LM{shank}.AllUnits ; fieldnames(getfield(LM{shank}.tableClusterData,ClusterTypes{i}))];
            ClusterUnitSpikes = cellfun(@(x) getfield(getfield(LM{shank}.tableClusterData,ClusterTypes{i}),x) , fieldnames(getfield(LM{shank}.tableClusterData,ClusterTypes{i})),'UniformOutput', false);
            LM{shank}.AllUnitSpikeTimes = [LM{shank}.AllUnitSpikeTimes, ClusterUnitSpikes'];
        end
        if strcmp(ClusterTypes{i},'good')
            LM{shank}.SingleUnits = [LM{shank}.SingleUnits ; fieldnames(getfield(LM{shank}.tableClusterData,ClusterTypes{i}))];
            ClusterUnitSpikes = cellfun(@(x) getfield(getfield(LM{shank}.tableClusterData,ClusterTypes{i}),x) , fieldnames(getfield(LM{shank}.tableClusterData,ClusterTypes{i})),'UniformOutput', false);
            LM{shank}.SingleUnitSpikeTimes = [LM{shank}.SingleUnitSpikeTimes, ClusterUnitSpikes'];
            
        end
    end
end
%%%V1
for shank = 1:length(V1)
    
    [V1{shank}.tableSpikeData, V1{shank}.tableClusterData, V1{shank}.strParams] = loadKiloSortedSpikes(V1{shank}.path2spikes);
    ClusterTypes = fieldnames(V1{shank}.tableClusterData);
    V1{shank}.AllUnits=[];
    V1{shank}.SingleUnits=[];
    V1{shank}.tableClusterData.all=[];
    V1{shank}.AllUnitSpikeTimes=cell(1,length(V1{shank}.AllUnits));
    V1{shank}.SingleUnitSpikeTimes=cell(1,length(V1{shank}.AllUnits));
    for i=1:length(ClusterTypes)
        if ~strcmp(ClusterTypes{i},'noise')
            V1{shank}.AllUnits = [V1{shank}.AllUnits ; fieldnames(getfield(V1{shank}.tableClusterData,ClusterTypes{i}))];
            ClusterUnitSpikes = cellfun(@(x) getfield(getfield(V1{shank}.tableClusterData,ClusterTypes{i}),x) , fieldnames(getfield(V1{shank}.tableClusterData,ClusterTypes{i})),'UniformOutput', false);
            V1{shank}.AllUnitSpikeTimes = [V1{shank}.AllUnitSpikeTimes, ClusterUnitSpikes'];
            
        end
        if strcmp(ClusterTypes{i},'good')
            V1{shank}.SingleUnits = [V1{shank}.SingleUnits ; fieldnames(getfield(V1{shank}.tableClusterData,ClusterTypes{i}))];
            ClusterUnitSpikes = cellfun(@(x) getfield(getfield(V1{shank}.tableClusterData,ClusterTypes{i}),x) , fieldnames(getfield(V1{shank}.tableClusterData,ClusterTypes{i})),'UniformOutput', false);
            V1{shank}.SingleUnitSpikeTimes = [V1{shank}.SingleUnitSpikeTimes, ClusterUnitSpikes'];
        end
    end
end

%% depths
%LM
for shank = 1:length(LM)
    temps=readNPY([LM{shank}.path2spikes,'templates.npy']);
    
    map=[LM{shank}.tableSpikeData.spike_templates LM{shank}.tableSpikeData.spike_clusters]; % this is not the best way. Just gets one correposnding cluster
    clusind = nan(1,max(LM{shank}.tableSpikeData.spike_clusters));
    clusind(map(:,2)+1) = map(:,1)+1;
    
    unitnumbers=cellfun(@(x) str2num(x(end-2:end)),LM{shank}.AllUnits);
    LM{shank}.AllUnitsDepths=arrayfun(@(x) giveDepth(temps,x), clusind(unitnumbers+1));
    
    if ~isempty(LM{shank}.SingleUnits)
        unitnumbers=cellfun(@(x) str2num(x(end-2:end)),LM{shank}.SingleUnits);
        LM{shank}.SingleUnitsDepths=arrayfun(@(x) giveDepth(temps,x), clusind(unitnumbers+1));
    end
end

%V1
for shank = 1:length(V1)
    temps=readNPY([V1{shank}.path2spikes,'templates.npy']);
    
     map=[V1{shank}.tableSpikeData.spike_templates V1{shank}.tableSpikeData.spike_clusters];
    clusind = nan(1,max(V1{shank}.tableSpikeData.spike_clusters));
    clusind(map(:,2)+1) = map(:,1)+1;
    
    unitnumbers=cellfun(@(x) str2num(x(end-2:end)),V1{shank}.AllUnits);
    V1{shank}.AllUnitsDepths=arrayfun(@(x) giveDepth(temps,x), clusind(unitnumbers+1));
    
    if ~isempty(V1{shank}.SingleUnits)
        unitnumbers=cellfun(@(x) str2num(x(end-2:end)),V1{shank}.SingleUnits);
        V1{shank}.SingleUnitsDepths=arrayfun(@(x) giveDepth(temps,x), clusind(unitnumbers+1));
    end
    
end
end
%%
function Depth = giveDepth(temps,x)
[~,Depth]=max(max(abs(squeeze(temps(x,:,:)))));
end
