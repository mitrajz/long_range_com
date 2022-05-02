
clear indss0_go indss0_nogo
edgestep = 20*30;
%min_n_trials_per_delay = 10;
onlynogrooming = 1;
poolalllags = 1;

if ~onlystableunits
    unitname = 'V1_Units.MUSpikes';
    units=(eval(unitname));
else
    units=([cell2mat(V1{1}.SingleUnitSpikeTimes(find(V1{1}.SingleUnitDriftMeasure_stable<driftthreshold))');...
        cell2mat(V1{2}.SingleUnitSpikeTimes(find(V1{2}.SingleUnitDriftMeasure_stable<driftthreshold))')]);
    
end
% 
% fstab = figure;
% units = ([cell2mat(V1{1}.SingleUnitSpikeTimes(9)')]);


rasterbinned=nan(size(PAllOn,1),(length(1:edgestep:size(PAllOn,2))-1));
for i=1:1:size(PAllOn,1)
    edges = PAllOn(i,1):edgestep:PAllOn(i,end);
    [rasterbinned(i,:),~] = histcounts(units,edges);
end
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

countspikesin = [floor(size(rasterbinned,2)/2) floor(size(rasterbinned,2)*3/4)];

figure(fstab);
ax=subplot(3,3,4+3*onlystableunits);plot(indss0_go,sum(rasterbinned(indss0_go,countspikesin(1):countspikesin(2)),2),'g.-')
hold on;plot(indss0_nogo,sum(rasterbinned(indss0_nogo,countspikesin(1):countspikesin(2)),2),'r.-');hold off
if ~onlystableunits
    ax.Title.String = 'MU spikes over trials';
else
    ax.Title.String = 'spikes of units that are stable over stable trials over trials';
end
figure(fstab);
subplot(3,3,5+3*onlystableunits);imagesc([],indss0_go,rasterbinned(indss0_go,:));
subplot(3,3,6+3*onlystableunits);imagesc([],indss0_nogo,rasterbinned(indss0_nogo,:))


% allspikesgo=sum(rasterbinned(indss0_go,countspikesin(1):countspikesin(2)),2);
% allspikesnogo=sum(rasterbinned(indss0_nogo,countspikesin(1):countspikesin(2)),2);
% unitvar1 = mean([std(allspikesgo)/mean(allspikesgo), std(allspikesnogo)/mean(allspikesnogo)])
% godrift = corrcoef(allspikesgo,indss0_go);
% nogodrift = corrcoef(allspikesnogo,indss0_nogo);
% unitvar2 = abs(mean([godrift(1,2) nogodrift(1,2)]))

%%
if false
    %rasterbinnedud = flipud(rasterbinned);
    % test=zeros(size(cormat));
    % indmat= sub2ind(size(cormat), 1:length(indss0_go)-1, 2:length(indss0_go));
    % test(indmat) = 1;
    % figure;imagesc(indmat)
    %%%%%%%%%%% cov or corcoeff?
    %%%%% correlation with 1 trial after
    cormat_go=(corrcoef(rasterbinned(indss0_go,:)'));
    indmat_go= sub2ind(size(cormat_go), 1:length(indss0_go)-1, 2:length(indss0_go));
    cormat_nogo=(corrcoef(rasterbinned(indss0_nogo,:)'));
    indmat_nogo= sub2ind(size(cormat_nogo), 1:length(indss0_nogo)-1, 2:length(indss0_nogo));
    figure;plot(indss0_go(1:end-1),cormat_go(indmat_go),'g.-');
    hold on;plot(indss0_nogo(1:end-1),cormat_nogo(indmat_nogo),'r.-');
    %%%%% correlation with average
    figure;plot(indss0_go(1:end),mean(cormat_go,1),'g.-');
    hold on;plot(indss0_nogo(1:end),mean(cormat_nogo,1),'r.-');
    avavcorcor = mean([mean(mean(cormat_nogo)) mean(mean(cormat_go))])
end
