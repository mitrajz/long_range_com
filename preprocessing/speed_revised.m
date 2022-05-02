% modifications from simple behavior:
% - trial by trial windowed, instead of the whole trace: less data, no need
% to include Rslag
% - binned instead of smoothed
% - position only goes forward (delta_position<0 is naned)
x=1;
speedWindowms =10;
thl=500;

lickch = n_ephys_channles + n_aux_channels + lick_channel;
RSch=n_ephys_channles + n_aux_channels + running_channel;


Pos = h5read(path2kwd{x},kwdrecording,[RSch 1],[1 inf]);
Rslag{x} = find(Pos > 0.99*max(Pos),1) - find(Pos < 1.1*min(Pos),1);

if Rslag{x}<15000 || Rslag{x}>45000
    Rslag{x} = 30000;
end

gain=h5read(path2kwd{x},'/recordings/0/application_data/channel_bit_volts');

%%% lick
alllick{x} = h5read(path2kwd{x},kwdrecording,[lickch 1],[1 inf]);
alllick{x}(alllick{x}<thl)=0;
alllick{x}(alllick{x}>thl)=1;

speedticks = nan(size(PAllOn{x},1),size(PAllOn{x},2)-1);
numspeedbins = (Window+WindowL)/(speedWindowms*(samplingf/1000));
speed= nan(size(PAllOn{x},1),numspeedbins);
licks = nan(size(PAllOn{x}));

for trialind=1:1:size(PAllOn{x},1)
 
    if ((PAllOn{x}(trialind,end) + Rslag{x}) - length(Pos)) <= 0
        %%%% Rs
        trialPos=double(Pos(PAllOn{x}(trialind,:) + Rslag{x})')*double(gain(RSch));
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
        %figure;plot(RSs{x}(PAllOn{x}(trialind,:) + Rslag{x})')
        %%%% licks
        licks(trialind,:) = alllick{x}(PAllOn{x}(trialind,:) + Rslag{x});
    end
end
%% grooming

nogroomingind=cell(size(expname));
nogroomingind{x} = [];
for trialind=1:1:size(PAllOn{x},1) 
    if ( (nanmean(speed(trialind,:))> grooming_speed_threshold) || (nanmean(licks(trialind,:)) < grooming_lick_threshold) )
        nogroomingind{x}(end+1) = trialind;
    end
end


%% plots

figure;beh_f1=gcf;beh_f1.Position = [400 400 1000 800];beh_f1.Name= 'no excluding grooming';
subplot(2,2,1);ax1=imagesc(((-Window+1)/samplingf:(speedWindowms/1000):Window/samplingf),[],...
    [speed(gotrialind{x},:);nan(1,size(speed,2));speed(nogotrialind{x},:)]);
ax1.Parent.XLabel.String='time(s)';
ax1.Parent.YLabel.String='trials';
ax1.Parent.Title.String=sprintf('Running speed (binned %d ms) aligned to stimulus onset',speedWindowms);


subplot(2,2,2);ax1=imagesc([-Window/samplingf,Window/samplingf],[],...
    [licks(gotrialind{x},:);nan(1,size(licks,2));licks(nogotrialind{x},:)]);
ax1.Parent.XLabel.String='time(s)';
ax1.Parent.YLabel.String='trials';
ax1.Parent.Title.String='licks aligned to stimulus onset';


ax1=subplot(2,2,3);plot(((-Window+1)/samplingf:(speedWindowms/1000):Window/samplingf),speed(gotrialind{x},:));
hold on;plot(((-Window+1)/samplingf:(speedWindowms/1000):Window/samplingf),nanmean(speed(gotrialind{x},:)),'k','LineWidth',2)

ax2=subplot(2,2,4);plot(((-Window+1)/samplingf:(speedWindowms/1000):Window/samplingf),speed(nogotrialind{x},:));
hold on;plot(((-Window+1)/samplingf:(speedWindowms/1000):Window/samplingf),nanmean(speed(nogotrialind{x},:)),'k','LineWidth',2)

ax1.XLabel.String='time(s)';
ax1.YLabel.String='running speed (cm/s)';
ax1.Title.String='Running speed go trials';
ax2.XLabel.String='time(s)';
ax2.YLabel.String='running speed (cm/s)';
ax2.Title.String='Running speed nogo trials';
ax1.YLim=[0 150];
ax2.YLim=[0 150];

%%%%% exclude grooming

figure;beh_f2=gcf;beh_f2.Position = [400 400 1000 800];beh_f2.Name= 'excluding grooming';
subplot(2,2,1);ax1=imagesc(((-Window+1)/samplingf:(speedWindowms/1000):Window/samplingf),[],...
    [speed(intersect(gotrialind{x},nogroomingind{x}),:);nan(1,size(speed,2));speed(intersect(nogotrialind{x},nogroomingind{x}),:)]);
ax1.Parent.XLabel.String='time(s)';
ax1.Parent.YLabel.String='trials';
ax1.Parent.Title.String=sprintf('Running speed (binned %d ms) aligned to stimulus onset',speedWindowms);



subplot(2,2,2);ax1=imagesc([-Window/samplingf,Window/samplingf],[],...
    [licks(intersect(gotrialind{x},nogroomingind{x}),:);nan(1,size(licks,2));licks(intersect(nogotrialind{x},nogroomingind{x}),:)]);
ax1.Parent.XLabel.String='time(s)';
ax1.Parent.YLabel.String='trials';
ax1.Parent.Title.String='licks aligned to stimulus onset';


ax1=subplot(2,2,3);plot(((-Window+1)/samplingf:(speedWindowms/1000):Window/samplingf),speed(intersect(gotrialind{x},nogroomingind{x}),:));
hold on;plot(((-Window+1)/samplingf:(speedWindowms/1000):Window/samplingf),nanmean(speed(intersect(gotrialind{x},nogroomingind{x}),:)),'k','LineWidth',2)

ax2=subplot(2,2,4);plot(((-Window+1)/samplingf:(speedWindowms/1000):Window/samplingf),speed(intersect(nogotrialind{x},nogroomingind{x}),:));
hold on;plot(((-Window+1)/samplingf:(speedWindowms/1000):Window/samplingf),nanmean(speed(intersect(nogotrialind{x},nogroomingind{x}),:)),'k','LineWidth',2)


ax1.XLabel.String='time(s)';
ax1.YLabel.String='running speed (cm/s)';
ax1.Title.String='Running speed go trials';
ax2.XLabel.String='time(s)';
ax2.YLabel.String='running speed (cm/s)';
ax2.Title.String='Running speed nogo trials';
ax1.YLim=[0 150];
ax2.YLim=[0 150];






