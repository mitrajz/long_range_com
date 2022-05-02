% targetcells
function targetcell = subselectdynamiccells(targetcell,subratio)
% after assigning targetcells and notargetcells 

% get average of all cells over time
cellagAv = nan(length(targetcell),8);
for celli = 1:length(targetcell)
    for lag = 1:8
       cellagAv(celli,lag) = nanmean(targetcell{celli}.laAbs.go{lag});
    end
    
end


%allcellagAv = nanmean(cellagAv,1); %size: 1*8 0 average of all cels per lag

ZscellagAv = nan(length(targetcell),8);
for lag = 1:8
       ZscellagAv(find(~isnan(cellagAv(:,lag))),lag) = ...
       zscore(cellagAv(find(~isnan(cellagAv(:,lag))),lag));
end

% cell indices ordered for highest variance over time (nans are in the end)
if subratio<0.5
    [~,ind] = sort(std(ZscellagAv'),'descend','MissingPlacement','last');
else
    [~,ind] = sort(std(ZscellagAv'),'ascend','MissingPlacement','last');
    subratio = 1 - subratio;
end

targetcell = targetcell(ind(1:round(subratio*numel(ind))));
