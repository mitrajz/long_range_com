

speedWindowms =10;
checkpos = 0;
speedticks = nan;
speed =nan;


valvech = n_ephys_channles + n_aux_channels + valve_channel;
lickch = n_ephys_channles + n_aux_channels + lick_channel;
RSch=n_ephys_channles + n_aux_channels + running_channel;

if checkpos
    Pos = h5read(path2kwd,kwdrecording,[RSch 1],[1 inf]);
    Rslag = find(Pos > 0.99*max(Pos),1) - find(Pos < 1.1*min(Pos),1);
    
    if Rslag<15000 || Rslag>45000
        Rslag = 30000;
    end
end

gain=h5read(path2kwd,'/recordings/0/application_data/channel_bit_volts');

%%% lick and valve
alllick = h5read(path2kwd,kwdrecording,[lickch 1],[1 inf]);
allvalve = h5read(path2kwd,kwdrecording,[valvech 1],[1 inf]);

%% normalize lick and valves
rectalllick = abs(alllick-prctile(alllick,50));

lickthreshold = prctile(rectalllick,lickthprctile);
valvethreshold = min(allvalve) + (max(allvalve)-min(allvalve))/2;%prctile(allvalve,99.9);%99.5

lick_th = rectalllick;
lick_th(rectalllick<lickthreshold) = 0;
lick_th(rectalllick>=lickthreshold) = 1;

valve_th = allvalve;
valve_th(allvalve<valvethreshold) = 0;
valve_th(allvalve>=valvethreshold) = 1;

%%%
if checkpos
    speedticks = nan(size(PAllOn,1),size(PAllOn,2)-1);
    numspeedbins = (Window+WindowL)/(speedWindowms*(samplingf/1000));
    speed= nan(size(PAllOn,1),numspeedbins);
end
licks = nan(size(PAllOn,1),floor(size(PAllOn,2)));
valves = nan(size(PAllOn,1),floor(size(PAllOn,2)));

for trialind=1:1:size(PAllOn,1)

    if checkpos
        if ((PAllOn(trialind,end) + Rslag) - length(Pos)) <= 0
            %%%% Rs
            trialPos=double(Pos(PAllOn(trialind,:) + Rslag)')*double(gain(RSch));
            %figure;plot(trialPos)
            speedticks(trialind,:)=double(diff(trialPos));
            speedticks(trialind,(speedticks(trialind,:))<0)=nan; % negative: going back or resetting position
            speedticks(trialind,(speedticks(trialind,:)<0.01))=0; % noise around a contant value
            speedticks(trialind,(speedticks(trialind,:)>0.5)) =nan; % >0.5cm steps
            
            for n=1:numspeedbins
                speed(trialind,n) = nansum(speedticks(trialind,(1+(n-1)*(speedWindowms*(samplingf/1000))):(n*(speedWindowms*(samplingf/1000)))))...
                    /(speedWindowms/1000);% in cm/s
            end
            %figure;plot(speed)
            %figure;plot(RSs(PAllOn(trialind,:) + Rslag)')
            %%%% licks
        end
    end
    licks(trialind,:) = lick_th(floor(PAllOn(trialind,:)));
    valves(trialind,:) = valve_th(floor(PAllOn(trialind,:)));
    
end

%% first licks and comparison to valve
lickcontforsample = lickcontforms*30;
firstlicksample = nan(1,size(PAllOn,1));
firstvalvesample = nan(1,size(PAllOn,1));
for trialind=1:1:size(PAllOn,1)
   
    fl=find(single(licks(trialind,floor(size(licks,2)/2):end)) == 1,1);
    fv=find(single(valves(trialind,floor(size(licks,2)/2):end)) == 1,1);

    if numel(fl) && fl>0
        firstlicksample(trialind) = fl;
    end
    
    if numel(fv) && fv>0
        firstvalvesample(trialind) = fv;
    end
end
f_lick1=figure;plot(firstlicksample/30,'.')
hold on;plot(firstvalvesample/30,'r.')
f_lick2=figure;histogram(firstlicksample/30, 100,'Normalization','count','FaceColor','b')
hold on;histogram(firstvalvesample/30, 100,'Normalization','count','FaceColor','r')
respwindowms = [100 650]; 
repwindowsample = 30*respwindowms;


correctnogotrialind = [];
correctgotrialind = [];
incorrectnogotrialind = [];
incorrectgotrialind = [];

FA=0;
Miss=0;

indss_nogo = nogotrialind;
indss_go = gotrialind;

for i = 1:length(indss_nogo)
    if (firstlicksample(indss_nogo(i))>repwindowsample(1)) & (firstlicksample(indss_nogo(i))<repwindowsample(2))
        FA=FA+1;
        incorrectnogotrialind(end+1) = indss_nogo(i);
    elseif (firstlicksample(indss_nogo(i))<repwindowsample(1))

    else
        correctnogotrialind(end+1) = indss_nogo(i);
    end
end

for i = 1:length(indss_go)
    if (firstlicksample(indss_go(i))>repwindowsample(1)) & (firstlicksample(indss_go(i))<repwindowsample(2))
        correctgotrialind(end+1) = indss_go(i);
    % if first lick before resp window, neither miss nor hit 
    elseif (firstlicksample(indss_go(i))<repwindowsample(1))
    else
        Miss=Miss+1;
        incorrectgotrialind(end+1) = indss_go(i);
    end

end

FA= 100*FA/length(indss_nogo)
Miss=100*Miss/length(indss_go)

correcttrials = zeros(1,length(size(PAllOn,1)));
incorrecttrials = zeros(1,length(size(PAllOn,1)));
correcttrials([correctnogotrialind correctgotrialind]) = 1;
incorrecttrials([incorrectnogotrialind incorrectgotrialind]) = 1;
f_behovertime = figure;plot(movmean(correcttrials, [0 20]));

%% grooming

nogroomingind=cell(size(expname));
nogroomingind = [];
for trialind=1:1:size(PAllOn,1) % grooming calculated only for prestimulus window
    if ((nanmean(licks(trialind,1:floor(size(licks,2)/2))) < grooming_lick_threshold) )
        nogroomingind(end+1) = trialind;
    end
end


%% plots

figure;
beh_f1=gcf;beh_f1.Name= 'including grooming';
beh_f1.Units = 'normalized';
beh_f1.Position = [0.3 0.4 0.5 0.4];
subplot(1,2,1);
ax1=imagesc(((-Window+1)/samplingf:(speedWindowms/1000):Window/samplingf),[],...
    [licks(gotrialind,:);nan(1,size(licks,2));licks(nogotrialind,:)]);
colormap(ax1.Parent,'gray')
ax1.Parent.XLabel.String='time(s)';
ax1.Parent.YLabel.String='trials';
ax1.Parent.Title.String='licks aligned to stimulus onset';
subplot(1,2,2);
ax1=imagesc(((-Window+1)/samplingf:(speedWindowms/1000):Window/samplingf),[],...
    [valves(gotrialind,:);nan(1,size(valves,2));valves(nogotrialind,:)],'AlphaData',0.5);
colormap(ax1.Parent,'gray')
ax1.Parent.XLabel.String='time(s)';
ax1.Parent.YLabel.String='trials';
ax1.Parent.Title.String='valve aligned to stimulus onset';

figure;
beh_f2=gcf;beh_f2.Name= 'excluding grooming';
beh_f2.Units = 'normalized';
beh_f2.Position = [0.3 0.4 0.5 0.4];
subplot(1,2,1);
ax1=imagesc(((-Window+1)/samplingf:(speedWindowms/1000):Window/samplingf),[],...
    [licks(intersect(gotrialind,nogroomingind),:);nan(1,size(licks,2));licks(intersect(nogotrialind,nogroomingind),:)]);
colormap(ax1.Parent,'gray')
ax1.Parent.XLabel.String='time(s)';
ax1.Parent.YLabel.String='trials';
ax1.Parent.Title.String='licks aligned to stimulus onset';
subplot(1,2,2);
ax1=imagesc(((-Window+1)/samplingf:(speedWindowms/1000):Window/samplingf),[],...
    [valves(intersect(gotrialind,nogroomingind),:);nan(1,size(valves,2));valves(intersect(nogotrialind,nogroomingind),:)],'AlphaData',0.5);
colormap(ax1.Parent,'gray')
ax1.Parent.XLabel.String='time(s)';
ax1.Parent.YLabel.String='trials';
ax1.Parent.Title.String='valve aligned to stimulus onset';
%% first lick vs. laserlag
% make a second inspection
min_n_trials_per_delay = 5;
cumsumstep=10;
onlynogrooming=0;
onlystableperiod=0;

normalizeforcorrect = 1;
figperfn = figure; figperfn.Name = 'normalized';
figperfn.Units = 'normalized';
figperfn.Position = [0.2 0.2 0.8 0.8];
makelaserlagperfplots;

normalizeforcorrect = 0;
figperf = figure; figperf.Name = 'not normalized';
figperf.Units = 'normalized';
figperf.Position = [0.2 0.2 0.8 0.8];
makelaserlagperfplots;
