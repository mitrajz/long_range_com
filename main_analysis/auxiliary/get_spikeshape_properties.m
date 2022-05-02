% get spike amplitude and width based on prctmax, and returns by updating
% fields of V1cells, LMcells
function [LMcells,V1cells,Behcells,params] = get_spikeshape_properties(prctmax,LMcells,V1cells,Behcells,params)

for i=1:length(V1cells)
    V1cells{i}.spikeamplitude = min(V1cells{i}.spiketemplate);
    V1cells{i}.spikewidth = ...
    (find(V1cells{i}.spiketemplate<prctmax*V1cells{i}.spikeamplitude,1,'last')) - ...
    (find(V1cells{i}.spiketemplate<prctmax*V1cells{i}.spikeamplitude,1,'first'));
end


for i=1:length(LMcells)
    LMcells{i}.spikeamplitude = min(LMcells{i}.spiketemplate);
    LMcells{i}.spikewidth = ...
    (find(LMcells{i}.spiketemplate<prctmax*LMcells{i}.spikeamplitude,1,'last')) - ...
    (find(LMcells{i}.spiketemplate<prctmax*LMcells{i}.spikeamplitude,1,'first'));
end





