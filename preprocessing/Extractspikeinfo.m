function [V1,LM] = Extractspikeinfo(V1,LM)


%%% LM
for shank = 1:2
    [V1{shank}.tableSpikeData, V1{shank}.tableClusterData, V1{shank}.strParams] = loadKiloSortedSpikes(V1{shank}.path2spikes);
    if isfield(V1{shank}.tableClusterData,'good')
        V1{shank}.GoodUnits=fieldnames(V1{shank}.tableClusterData.good);
        V1{shank}.GUnitSpikeTimes=cell(1,length(V1{shank}.GoodUnits));
        for i=1:1:length(V1{shank}.GoodUnits)
            V1{shank}.GUnitSpikeTimes{i}=getfield(V1{shank}.tableClusterData.good,V1{shank}.GoodUnits{i}); % in oe samples
        end
    end
    % getting all excep for noise
    ClusterTypes = fieldnames(V1{shank}.tableClusterData);
    V1{shank}.AllUnits=[];
    V1{shank}.tableClusterData.all=[];
    V1{shank}.AllUnitSpikeTimes=cell(1,length(V1{shank}.AllUnits));
    for i=1:length(ClusterTypes)
        if ~strcmp(ClusterTypes{i},'noise')
            V1{shank}.AllUnits = [V1{shank}.AllUnits ; fieldnames(getfield(V1{shank}.tableClusterData,ClusterTypes{i}))];
            %V1{shank}.tableClusterData.all = [V1{shank}.tableClusterData.all ; getfield(V1{shank}.tableClusterData,ClusterTypes{i})];
            ClusterUnitSpikes = cellfun(@(x) getfield(getfield(V1{shank}.tableClusterData,ClusterTypes{i}),x) , fieldnames(getfield(V1{shank}.tableClusterData,ClusterTypes{i})),'UniformOutput', false);
            V1{shank}.AllUnitSpikeTimes = [V1{shank}.AllUnitSpikeTimes, ClusterUnitSpikes'];
        end
    end
end

%%% LM
for shank = 1:2
    [LM{shank}.tableSpikeData, LM{shank}.tableClusterData, LM{shank}.strParams] = loadKiloSortedSpikes(LM{shank}.path2spikes);
    if isfield(LM{shank}.tableClusterData,'good')
        LM{shank}.GoodUnits=fieldnames(LM{shank}.tableClusterData.good);
        LM{shank}.GUnitSpikeTimes=cell(1,length(LM{shank}.GoodUnits));
        for i=1:1:length(LM{shank}.GoodUnits)
            LM{shank}.GUnitSpikeTimes{i}=getfield(LM{shank}.tableClusterData.good,LM{shank}.GoodUnits{i}); % in oe samples
        end
    end
    % getting all excep for noise
    ClusterTypes = fieldnames(LM{shank}.tableClusterData);
    LM{shank}.AllUnits=[];
    LM{shank}.tableClusterData.all=[];
    LM{shank}.AllUnitSpikeTimes=cell(1,length(LM{shank}.AllUnits));
    for i=1:length(ClusterTypes)
        if ~strcmp(ClusterTypes{i},'noise')
            LM{shank}.AllUnits = [LM{shank}.AllUnits ; fieldnames(getfield(LM{shank}.tableClusterData,ClusterTypes{i}))];
            %V1{shank}.tableClusterData.all = [V1{shank}.tableClusterData.all ; getfield(V1{shank}.tableClusterData,ClusterTypes{i})];
            ClusterUnitSpikes = cellfun(@(x) getfield(getfield(LM{shank}.tableClusterData,ClusterTypes{i}),x) , fieldnames(getfield(LM{shank}.tableClusterData,ClusterTypes{i})),'UniformOutput', false);
            LM{shank}.AllUnitSpikeTimes = [LM{shank}.AllUnitSpikeTimes, ClusterUnitSpikes'];
        end
    end
end
