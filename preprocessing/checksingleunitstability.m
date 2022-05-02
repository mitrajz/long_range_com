%% number of spikes during stimulus is a better measure of single cell variability compared to correlations
% this version stability of units only takes into account the drift with time, 
% that is if there is a linear increase/decrease in the overal number of
% spikes during the stimulus with time. 

clear indss0_go indss0_nogo

edgestep = 20*30;
% min_n_trials_per_delay = 10;
onlynogrooming = 1;
poolalllags = 1;
plotsanitycheck = 0;


if onlynogrooming
    if poolalllags
        indss0_nogo=intersect(nogotrialind,nogroomingind);
        indss0_go=intersect(gotrialind,nogroomingind);
    else
        indss0_nogo=intersect(intersect(find(isnan(LaserDelay)),nogotrialind),nogroomingind);
        indss0_go=intersect(intersect(find(isnan(LaserDelay)),gotrialind),nogroomingind);
    end
else
    if poolalllags
        indss0_nogo=nogotrialind;
        indss0_go=gotrialind;
    else
        indss0_nogo= intersect(find(isnan(LaserDelay)),nogotrialind);
        indss0_go= intersect(find(isnan(LaserDelay)),gotrialind);
    end
end

if onlystableperiod
    indss0_go = intersect(indss0_go, stabletrialind);
    indss0_nogo = intersect(indss0_nogo, stabletrialind);
end


areanames = {'V1','LM'};
for areai = 1:2
    area = areanames{areai};
    for shank = 1:2
        TempElectrode = eval(sprintf('%s{%d}',area,shank));
        for unitnum = 1:length(TempElectrode.SingleUnitSpikeTimes)
            unit = TempElectrode.SingleUnitSpikeTimes{unitnum};
            
            rasterbinned=nan(size(PAllOn,1),(length(1:edgestep:size(PAllOn,2))-1));
            for i=1:1:size(PAllOn,1)
                edges = PAllOn(i,1):edgestep:PAllOn(i,end);
                [rasterbinned(i,:),~] = histcounts(unit,edges);
            end
            countspikesin = [floor(size(rasterbinned,2)/2) floor(size(rasterbinned,2)*3/4)];
            
            allspikesgo=sum(rasterbinned(indss0_go,countspikesin(1):countspikesin(2)),2);
            allspikesnogo=sum(rasterbinned(indss0_nogo,countspikesin(1):countspikesin(2)),2);
            % unitvar1 = mean([std(allspikesgo)/mean(allspikesgo), std(allspikesnogo)/mean(allspikesnogo)])
            godrift = corrcoef(allspikesgo,indss0_go);
            nogodrift = corrcoef(allspikesnogo,indss0_nogo);
            unitvar = abs(mean([godrift(1,2) nogodrift(1,2)]));
            if onlystableperiod
                eval(sprintf('%s{%d}.SingleUnitDriftMeasure_stable(unitnum) = unitvar;',area,shank));
            elseif ~onlystableperiod
                eval(sprintf('%s{%d}.SingleUnitDriftMeasure(unitnum) = unitvar;',area,shank));
            end
            if plotsanitycheck               
                figure;subplot(2,2,1);imagesc(rasterbinned(indss0_go,:));subplot(2,2,2);imagesc(rasterbinned(indss0_nogo,:))
                subplot(2,2,[3 4]);plot(indss0_go,sum(rasterbinned(indss0_go,countspikesin(1):countspikesin(2)),2),'g.-')
                hold on;plot(indss0_nogo,sum(rasterbinned(indss0_nogo,countspikesin(1):countspikesin(2)),2),'r.-')
                title(sprintf('single unit drift measure:%f',unitvar));
                
            end
        end
    end
end






