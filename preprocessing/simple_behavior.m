% this is a simple visual-only behavior analysis

lickch = n_ephys_channles + n_aux_channels + lick_channel;
RSch=n_ephys_channles + n_aux_channels + running_channel;

Rslag = cell(size(expname));
lick = cell(size(expname));
RS = cell(size(expname));
RSs = cell(size(expname));

nogroomingind=cell(size(expname));

for x=1:1:length(expname)
    %%% RS    
    Pos = h5read(path2kwd{x},kwdrecording,[RSch 1],[1 inf]);
    RsstartRec = find(Pos > 0.9*max(Pos),1);
    Rslag{x} = find(Pos > 0.99*max(Pos),1) - find(Pos < 1.1*min(Pos),1);
    
    gain=h5read(path2kwd{x},'/recordings/0/application_data/channel_bit_volts');
    
    RS{x}=double(diff(Pos))*double(gain(RSch))/double((1/30000)); % should be cm per sec, no?
    RS{x}(RS{x}>10^8*double(gain(RSch))) =nan; % probably not corect:5
    RS{x}(RS{x}<0)=nan;
    RSs{x} = nanfastsmooth(RS{x},3001,1,1); %untrustworthy 1000
    
    %%% lick
    lick{x} = h5read(path2kwd{x},kwdrecording,[lickch 1],[1 inf]);
    thl=500;%1000%(max(lick{x})+min(lick{x}))/10;
    lick{x}(lick{x}<thl)=0;
    lick{x}(lick{x}>thl)=1;
    
end
%% grooming


for x=1:1:length(expname)
    

    nogroomingind{x}=find(((mean(RSs{x}(PAllOn{x}),2) < grooming_RSc_threshold) ...
        & (mean(lick{x}(PAllOn{x}),2)> grooming_lick_threshold)) == 0);
    
end


%% plots
figure;beh_f1=gcf;beh_f1.Position = [400 400 1000 800];beh_f1.Name= 'no excluding grooming';
subplot(2,2,1);ax1=imagesc([-Window/samplingf,Window/samplingf],[],[RSs{x}(PAllOn{x}(gotrialind{x}(1:end-1),:) + Rslag{x})',nan(size(PAllOn{x},2),1),RSs{x}(PAllOn{x}(nogotrialind{x}(1:end-1),:) + Rslag{x})']');
ax1.Parent.XLabel.String='time(s)';
ax1.Parent.YLabel.String='trials';
ax1.Parent.Title.String='Running speed aligned to stimulus onset';


subplot(2,2,2);ax1=imagesc([-Window/samplingf,Window/samplingf],[],[5*lick{x}(PAllOn{x}(gotrialind{x}(1:end-1),:) + Rslag{x})',ones(size(PAllOn{x},2),1),5*lick{x}(PAllOn{x}(nogotrialind{x}(1:end-1),:) + Rslag{x})']');
ax1.Parent.XLabel.String='time(s)';
ax1.Parent.YLabel.String='trials';
ax1.Parent.Title.String='licks aligned to stimulus onset';


ax1=subplot(2,2,3);plot((-WindowL:Window)/samplingf,RSs{x}(PAllOn{x}(gotrialind{x}(1:end-1),:) + Rslag{x})');
ax2=subplot(2,2,4);plot((-WindowL:Window)/samplingf,RSs{x}(PAllOn{x}(nogotrialind{x}(1:end-1),:) + Rslag{x})');
ax1.XLabel.String='time(s)';
ax1.YLabel.String='running speed (cm/s)';
ax1.Title.String='Running speed go trials';
ax2.XLabel.String='time(s)';
ax2.YLabel.String='running speed (cm/s)';
ax2.Title.String='Running speed nogo trials';

%%%%% exclude grooming

figure;beh_f2=gcf;beh_f2.Position = [400 400 1000 800];beh_f2.Name= 'excluding grooming';
subplot(2,2,1);ax1=imagesc([-Window/samplingf,Window/samplingf],[],[RSs{x}(PAllOn{x}(intersect(gotrialind{x}(1:end-1),nogroomingind{x}),:) + Rslag{x})',nan(size(PAllOn{x},2),1),RSs{x}(PAllOn{x}(intersect(nogotrialind{x}(1:end-1),nogroomingind{x}),:) + Rslag{x})']');
ax1.Parent.XLabel.String='time(s)';
ax1.Parent.YLabel.String='trials';
ax1.Parent.Title.String='Running speed aligned to stimulus onset';


subplot(2,2,2);ax1=imagesc([-Window/samplingf,Window/samplingf],[],[5*lick{x}(PAllOn{x}(intersect(gotrialind{x}(1:end-1),nogroomingind{x}),:) + Rslag{x})',ones(size(PAllOn{x},2),1),5*lick{x}(PAllOn{x}(intersect(nogotrialind{x}(1:end-1),nogroomingind{x}),:) + Rslag{x})']');
ax1.Parent.XLabel.String='time(s)';
ax1.Parent.YLabel.String='trials';
ax1.Parent.Title.String='licks aligned to stimulus onset';


ax1=subplot(2,2,3);plot((-WindowL:Window)/samplingf,RSs{x}(PAllOn{x}(intersect(gotrialind{x}(1:end-1),nogroomingind{x}),:) + Rslag{x})');
ax2=subplot(2,2,4);plot((-WindowL:Window)/samplingf,RSs{x}(PAllOn{x}(intersect(nogotrialind{x}(1:end-1),nogroomingind{x}),:) + Rslag{x})');
ax1.XLabel.String='time(s)';
ax1.YLabel.String='running speed (cm/s)';
ax1.Title.String='Running speed go trials';
ax2.XLabel.String='time(s)';
ax2.YLabel.String='running speed (cm/s)';
ax2.Title.String='Running speed nogo trials';







