function [V1cells,LMcells] = SubselectExpType(V1cells,LMcells,exp2analyse)

% choosing V1 cells from only FB or FF experiments
V1cells = V1cells(find(cellfun(@(x) strcmp(x.exptype,exp2analyse),V1cells)));
% choosing LM cells from only FB or FF experiments
LMcells = LMcells(find(cellfun(@(x) strcmp(x.exptype,exp2analyse),LMcells)));
%%%%%%%%% go and nogo laser delay distributions across all animals
figure;histogram(reshape(cell2mat(transpose(cellfun(@(x) x.smb_centers.go,LMcells,'UniformOutput',0))),1,[]),200,'FaceColor','g','FaceAlpha',0.3);
hold on;histogram(reshape(cell2mat(transpose(cellfun(@(x) x.smb_centers.nogo,LMcells,'UniformOutput',0))),1,[]),200,'FaceColor','r','FaceAlpha',0.3);

figure;histogram(reshape(cell2mat(transpose(cellfun(@(x) x.smb_centers.go,V1cells,'UniformOutput',0))),1,[]),200,'FaceColor','g','FaceAlpha',0.3);
hold on;histogram(reshape(cell2mat(transpose(cellfun(@(x) x.smb_centers.nogo,V1cells,'UniformOutput',0))),1,[]),200,'FaceColor','r','FaceAlpha',0.3);

