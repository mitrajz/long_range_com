function [allhitrate,allmissrate,allfarate,allcrrate,allxpoints] = ...
    behsummaryhelper_rates(gotrialind,nogotrialind,...
    correctgotrialind,incorrectgotrialind,correctnogotrialind,incorrectnogotrialind,...
    stabletrialind,nogroomingind,LaserDelayBinned,...
    min_n_trials_per_delay,onlynogrooming,onlystableperiod,centers)



gotrialind_w = gotrialind;
nogotrialind_w = nogotrialind;
correctgotrialind_w = correctgotrialind;
incorrectgotrialind_w = incorrectgotrialind;
correctnogotrialind_w = correctnogotrialind;
incorrectnogotrialind_w = incorrectnogotrialind;


if onlystableperiod
    gotrialind_w = intersect(gotrialind_w,stabletrialind);
    nogotrialind_w = intersect(nogotrialind_w,stabletrialind);
    correctgotrialind_w = intersect(correctgotrialind_w,stabletrialind);
    incorrectgotrialind_w = intersect(incorrectgotrialind_w,stabletrialind);
    correctnogotrialind_w = intersect(correctnogotrialind_w,stabletrialind);
    incorrectnogotrialind_w = intersect(incorrectnogotrialind_w,stabletrialind);
end
if onlynogrooming
    gotrialind_w = intersect(gotrialind_w,nogroomingind);
    nogotrialind_w = intersect(nogotrialind_w,nogroomingind);
    correctgotrialind_w = intersect(correctgotrialind_w,nogroomingind);
    incorrectgotrialind_w = intersect(incorrectgotrialind_w,nogroomingind);
    correctnogotrialind_w = intersect(correctnogotrialind_w,nogroomingind);
    incorrectnogotrialind_w = intersect(incorrectnogotrialind_w,nogroomingind);
end


allhitrate = [];
allmissrate = [];
allxpoints = [];
allcrrate = [];
allfarate = [];
%%%%%%%%%%%%%%% baseline 
% go
indss0 = intersect(find(isnan(LaserDelayBinned)),gotrialind_w);

allhitrate(end+1) = length(intersect(correctgotrialind_w,indss0))/...
    (length(intersect(correctgotrialind_w,indss0))+length(intersect(incorrectgotrialind_w,indss0)));
allmissrate(end+1) = length(intersect(incorrectgotrialind_w,indss0))/...
    (length(intersect(correctgotrialind_w,indss0))+length(intersect(incorrectgotrialind_w,indss0)));
% nogo
indss0= intersect(find(isnan(LaserDelayBinned)),nogotrialind_w);
allcrrate(end+1) = length(intersect(correctnogotrialind_w,indss0))/...
    (length(intersect(correctnogotrialind_w,indss0))+length(intersect(incorrectnogotrialind_w,indss0)));
allfarate(end+1) = length(intersect(incorrectnogotrialind_w,indss0))/...
    (length(intersect(correctnogotrialind_w,indss0))+length(intersect(incorrectnogotrialind_w,indss0)));

%%%%%%%%%%%%%%% laser
for ldelay =1:length(centers)
    
    if length(find(LaserDelayBinned == (ldelay))) > min_n_trials_per_delay
        % go
        indss= intersect(find(LaserDelayBinned == (ldelay)),gotrialind_w);
        allhitrate(end+1) = length(intersect(correctgotrialind_w,indss))/...
            (length(intersect(correctgotrialind_w,indss))+length(intersect(incorrectgotrialind_w,indss)));
        allmissrate(end+1) = length(intersect(incorrectgotrialind_w,indss))/...
            (length(intersect(correctgotrialind_w,indss))+length(intersect(incorrectgotrialind_w,indss)));
        % nogo
        indss= intersect(find(LaserDelayBinned == (ldelay)),nogotrialind_w); 
        allcrrate(end+1) = length(intersect(correctnogotrialind_w,indss))/...
            (length(intersect(correctnogotrialind_w,indss))+length(intersect(incorrectnogotrialind_w,indss)));
        allfarate(end+1) = length(intersect(incorrectnogotrialind_w,indss))/...
            (length(intersect(correctnogotrialind_w,indss))+length(intersect(incorrectnogotrialind_w,indss)));
        
        % xpoints
        allxpoints(end+1) = centers(ldelay);
        
    end
end







