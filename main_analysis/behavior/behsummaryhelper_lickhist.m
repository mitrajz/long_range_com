function   [lickcumsum,licklatencydist,licklatencydist_fa] = behsummaryhelper_lickhist(gotrialind,nogotrialind,...
    correctgotrialind,incorrectgotrialind,correctnogotrialind,incorrectnogotrialind,...
    stabletrialind,nogroomingind,LaserDelayBinned,...
    min_n_trials_per_delay,onlynogrooming,onlystableperiod,...
    cumsumstep,normalizeforcorrect,respwindowms,put2nanoutsidewin,centers,firstlicksample,onlycorrect)


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

lickcumsum = [];
licklatencydist=cell(0,0);
licklatencydist_fa=cell(0,0);
%%%%%%%%%%%%%%% baseline
% go
if onlycorrect
    indss0 = intersect(intersect(find(isnan(LaserDelayBinned)),gotrialind_w),correctgotrialind_w);
else
    indss0 = intersect(find(isnan(LaserDelayBinned)),gotrialind_w);
end
indssfirstlicks=firstlicksample(indss0)/30;
if put2nanoutsidewin
    indssfirstlicks(find(indssfirstlicks<respwindowms(1))) =nan;
    indssfirstlicks(find(indssfirstlicks>=respwindowms(2))) =nan;
end
if normalizeforcorrect
    indssfirstlicks(find(isnan(indssfirstlicks))) = [];
end
if onlycorrect
    licklatencydist{end+1} = indssfirstlicks;
else
    licklatencydist{end+1} = indssfirstlicks(find(indssfirstlicks>respwindowms(1) & indssfirstlicks<respwindowms(2)));
end
sig=cumsum(histcounts(indssfirstlicks,respwindowms(1):cumsumstep:respwindowms(2),'Normalization','count'))/numel(indssfirstlicks);
if isnan(sum(sig)) % if no licks found in the window
    sig(find(isnan(sig))) = 0;
end
lickcumsum(end+1,1:length(sig)) = sig;

% nogo: 

indss0 = intersect(find(isnan(LaserDelayBinned)),incorrectnogotrialind_w);

indssfirstlicks=firstlicksample(indss0)/30;
if put2nanoutsidewin
    indssfirstlicks(find(indssfirstlicks<respwindowms(1))) =nan;
    indssfirstlicks(find(indssfirstlicks>=respwindowms(2))) =nan;
end
if normalizeforcorrect
    indssfirstlicks(find(isnan(indssfirstlicks))) = [];
end
licklatencydist_fa{end+1} = indssfirstlicks;


%%%%%%%%%%%%%%% laser
for ldelay =1:length(centers)
    
    if length(find(LaserDelayBinned == (ldelay))) > min_n_trials_per_delay
        % go
        if onlycorrect
            indss= intersect(intersect(find(LaserDelayBinned == (ldelay)),gotrialind_w),correctgotrialind_w);  
        else
            indss= intersect(find(LaserDelayBinned == (ldelay)),gotrialind_w);  
        end
        indssfirstlicks=firstlicksample(indss)/30;
        if put2nanoutsidewin
            indssfirstlicks(find(indssfirstlicks<respwindowms(1))) =nan;
            indssfirstlicks(find(indssfirstlicks>=respwindowms(2))) =nan;
        end
        if normalizeforcorrect
            indssfirstlicks(find(isnan(indssfirstlicks))) = [];
        end
        if onlycorrect
            licklatencydist{end+1} = indssfirstlicks;
        else
            licklatencydist{end+1} = indssfirstlicks(find(indssfirstlicks>respwindowms(1) & indssfirstlicks<respwindowms(2)));
        end
        sig=cumsum(histcounts(indssfirstlicks,respwindowms(1):cumsumstep:respwindowms(2),'Normalization','count'))/numel(indssfirstlicks);
        lickcumsum(end+1,1:length(sig)) =sig;
        % nogo:
        
        indss= intersect(find(LaserDelayBinned == (ldelay)),incorrectnogotrialind_w);  
        
        indssfirstlicks=firstlicksample(indss)/30;
        if put2nanoutsidewin
            indssfirstlicks(find(indssfirstlicks<respwindowms(1))) =nan;
            indssfirstlicks(find(indssfirstlicks>=respwindowms(2))) =nan;
        end
        if normalizeforcorrect
            indssfirstlicks(find(isnan(indssfirstlicks))) = [];
        end     
        licklatencydist_fa{end+1} = indssfirstlicks;
    end
end