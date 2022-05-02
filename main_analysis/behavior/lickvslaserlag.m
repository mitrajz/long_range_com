% look at licks with respect to laser delay


%figure

% go trials:
laserlickgo = [];
laserlicknogo = [];
laserlickgotimes = [];
laserlicknogotimes = [];
% 
% hitrate = nan(9,8);
% crrate =  nan(9,8);
animaln =9;

ms1 = 50; % ms after visual onset
ms2 = 850;
firstpossiblefirstlick = floor(size(Behcells{animaln}.licksalltrials.go,2)/2) + ms1 /(edgestep/30);

Behcells{animaln}.firstlicksalltrials.go = zeros(size(Behcells{animaln}.licksalltrials.go));
Behcells{animaln}.firstlicksalltrials.nogo = zeros(size(Behcells{animaln}.licksalltrials.nogo));
for i=1:size(Behcells{animaln}.licksalltrials.go,1)
    Behcells{animaln}.firstlicksalltrials.go(i,firstpossiblefirstlick+find(Behcells{animaln}.licksalltrials.go(i,firstpossiblefirstlick:end)>0,1))=1;
    Behcells{animaln}.firstlicksalltrials.nogo(i,firstpossiblefirstlick+find(Behcells{animaln}.licksalltrials.nogo(i,firstpossiblefirstlick:end)>0,1))=1;
end

% define respwindow
respwindow = [firstpossiblefirstlick,  floor(size(Behcells{animaln}.licksalltrials.go,2)/2) + ms2 /(edgestep/30)];

for laserlag = 1:9
    % for 1 animal
    laserlickgo(end+1,:)=nanmean(Behcells{animaln}.firstlicksalltrials.go(Behcells{animaln}.ldelayindss.go{laserlag},:),1);
    laserlicknogo(end+1,:)=nanmean(Behcells{animaln}.firstlicksalltrials.nogo(Behcells{animaln}.ldelayindss.nogo{laserlag},:),1);
    
    licksinwindow = (sum(Behcells{animaln}.firstlicksalltrials.go(Behcells{animaln}.ldelayindss.go{laserlag},respwindow(1):respwindow(2)),2));
    hitrate(animaln,laserlag) = length(find(licksinwindow>0))/length(licksinwindow);
    
    licksinwindow = (sum(Behcells{animaln}.firstlicksalltrials.nogo(Behcells{animaln}.ldelayindss.nogo{laserlag},respwindow(1):respwindow(2)),2));
    crrate(animaln,laserlag) = 1-length(find(licksinwindow>0))/length(licksinwindow);
    
    
     laserlickgotimes(end+1,:)=[mean(find(laserlickgo(laserlag,:)>0)) std(find(laserlickgo(laserlag,:)>0))/sqrt(length(find(laserlickgo(laserlag,:)>0)))];
     laserlicknogotimes(end+1,:)=[mean(find(laserlicknogo(laserlag,:)>0)) std(find(laserlicknogo(laserlag,:)>0))/sqrt(length(find(laserlicknogo(laserlag,:)>0)))];

end

%% first index is indss (no laser)

figure;plot(hitrate(1:4,:)','*-')
xticklabels(arrayfun(@(x) sprintf('%d',x),[nan floor(V1cells{1}.smb_centers.go)],'UniformOutput',0))
figure;plot(crrate(1:4,:)','*-')
xticklabels(arrayfun(@(x) sprintf('%d',x),[nan floor(V1cells{1}.smb_centers.go)],'UniformOutput',0))