
%%
addpath(genpath('/mnt/data/Mitra/cache/codes/matlabkwikreader/'))
addpath(genpath('/mnt/data/Mitra/cache/codes/'))
%%
coef = -10
savefigs=0;
edgestep =40*30;%5
edgestepsmall = 1*30;
x=1;
min_n_trials_per_delay = 22;
onlygoandnogrooming = 1;
ask_for_save = 0;
unitname = 'V1{1}.SingleUnitSpikeTimes{5}'; 
units=(eval(unitname));


rasterbinned=nan(size(PAllOn,1),(length(1:edgestep:size(PAllOn,2))-1));
rasterbinnedsmall=nan(size(PAllOn,1),(length(1:edgestepsmall:size(PAllOn,2))-1));

for i=1:1:size(PAllOn,1)
    edges = PAllOn(i,1):edgestep:PAllOn(i,end);
    edgessmall = PAllOn(i,1):edgestepsmall:PAllOn(i,end);
    
    [rasterbinned(i,:),~] = histcounts(units,edges);
    [rasterbinnedsmall(i,:),~] = histcounts(units,edgessmall);
end
if onlygoandnogrooming
    indss0=intersect(intersect(find(isnan(LaserDelay)),gotrialind),nogroomingind);
else
    indss0= intersect(find(isnan(LaserDelay)),gotrialind);
end
trace=cell(0,0);
stemb = cell(0,0);
trace{end+1} = mean(rasterbinned(indss0,:),1);
stdmb0 = std(rasterbinned(indss0,:),0,1);
stemb{end+1} = stdmb0/sqrt(length(indss0));


count=1;
rastermat = rasterbinnedsmall(indss0,:);
x=cell(0,0);
y=cell(0,0);
l = cell(0,0);

for ldelay =1:length(centers)    
    if onlygoandnogrooming
        indss=intersect(intersect(find(LaserDelayBinned == (ldelay)),gotrialind),nogroomingind);

    else
        indss= intersect(find(LaserDelayBinned == (ldelay)),gotrialind);
    end
    if length(indss) > min_n_trials_per_delay
    
        trace{end+1} = mean(rasterbinned(indss,:),1);
        l{end+1} = length(indss);
        stdmb = std(rasterbinned(indss,:),0,1);
        stemb{end+1} = stdmb/sqrt(length(indss));
        
        
        x{end+1} = [centers(ldelay) centers(ldelay) centers(ldelay)+150 centers(ldelay)+150];
        y{end+1} = [size(rastermat,1)+1 size(rastermat,1)+length(indss)+2 ...
            size(rastermat,1)+length(indss)+2 size(rastermat,1)+1];
        rastermat = [rastermat;10*ones(1,size(rasterbinnedsmall,2));rasterbinnedsmall(indss,:)];
        count = count+1;
    end
end
f=figure;
f.Name = [animalname,'_',unitname];
f.Units = 'normalized';
f.Position = [0.3 0.3 0.5 0.8]


subplot(1,2,1);
timeoffset = .5*edgestepsmall/(samplingf/1000);
imagesc(timeoffset +((-Window):1*edgestepsmall:(Window-1))/(samplingf/1000),[],rastermat);
greymap=colormap('gray');
colormap(flipud(greymap));
set(gca,'CLim',[0 0.5]);
xlim([-100 600])

fig = gca;
hold on;
patch([0 0 500 500],[0 fig.YLim(2) fig.YLim(2) 0],[0 0 0],'FaceAlpha',0.05,'EdgeColor','none')

for i=1:length(x)

hold on;
patch(x{i},y{i},[0 0 1],'FaceAlpha',0.2,'EdgeColor','none')
end
box('off')

timeoffset = .5*edgestep/(samplingf/1000);
subplot(1,2,2)
for i=1:length(x)

hold on;
patch(x{i},[y{i}(1) y{i}(2) y{i}(3) y{i}(4)],[0 0 1],'FaceAlpha',0.05,'EdgeColor','none')
hold on;
%plot(((-Window):1*edgestep:(Window-1))/(samplingf/1000),-10*trace{i}+y{i}(1),'Color',[0 0.8 1],'LineWidth',1)
if i~= 1
shadedErrorBar(timeoffset + ((-Window):1*edgestep:(Window-1))/(samplingf/1000),...
    coef*trace{i}+y{i}(1),...
    coef*1*stemb{i},{'Color',[0 0 1],'Linewidth',1.5,'MarkerEdgeColor','none'},1)
end
    hold on;
%plot(((-Window):1*edgestep:(Window-1))/(samplingf/1000),-10*trace{1}+y{i}(1),'k','LineWidth',1)
shadedErrorBar(timeoffset + ((-Window):1*edgestep:(Window-1))/(samplingf/1000),...
    coef*trace{1}+y{i}(1),...
    coef*1*stemb{1},{'Color',[0.2 0.2 0.2],'Linewidth',1.5,'MarkerEdgeColor','none'},1)
if i==1
    hold on;
plot(timeoffset +((-Window):1*edgestep:(Window-1))/(samplingf/1000),coef*trace{i}+y{i}(1),'k','LineWidth',1)
end
xlim([-100 600])
end
hold on;
%plot(((-Window):1*edgestep:(Window-1))/(samplingf/1000),-10*trace{i+1}+y{i}(2),'Color',[0 0.8 1],'LineWidth',1)
shadedErrorBar(timeoffset + ((-Window):1*edgestep:(Window-1))/(samplingf/1000),...
    coef*trace{i+1}+y{i}(2),...
    coef*1*stemb{i+1},{'Color',[0 0 1],'Linewidth',1.5,'MarkerEdgeColor','none'},1)
hold on;
%plot(((-Window):1*edgestep:(Window-1))/(samplingf/1000),-10*trace{1}+y{i}(2),'k','LineWidth',1)
shadedErrorBar(timeoffset +((-Window):1*edgestep:(Window-1))/(samplingf/1000),...
    coef*trace{1}+y{i}(2),...
    coef*1*stemb{1},{'Color',[0.2 0.2 0.2],'Linewidth',1.5,'MarkerEdgeColor','none'},1)
xlim([-100 600])
ax = gca;
ax.YDir = 'reverse';
ax.YLim = fig.YLim;


hold on;
patch([0 0 500 500],[0 fig.YLim(2) fig.YLim(2) 0],[0 0 0],'FaceAlpha',0.05,'EdgeColor','none')
ax.YTick={};
box('off')



if savefigs
set(gcf,'Color','w');
saveas(f,[expname{1}(2:end),'_',f.Name,'.fig']);
export_fig(f,sprintf('%s',[f.Name,'.png']))
end
%%
cellnum=35;
figure;plot(V1cells{cellnum}.smb_centers.go,V1cells{cellnum}.smb.go,'g.-')
box('off')
ylabel('% Influence')
xlabel('Time(ms)')
ylim([-70 70])

randsmb80 = cell(1,8);
for nlag=1:8
    randrep=1000;
    randsmb80{nlag} = nan(1,randrep);
    for i=1:randrep

        pool=[V1cells{cellnum}.laAls.go{nlag};V1cells{cellnum}.laAbs.go{nlag}];
        pool=pool(randperm(length(pool)));
        nlasertrials= size(V1cells{cellnum}.laAls.go{nlag},1);
       
        randsmb80{nlag}(i) =  100*(mean(pool(1:nlasertrials)) - mean(pool(nlasertrials+1:end)))/...
             mean(pool(nlasertrials+1:end));
        
    end
end


hold on;patch('XData',[V1cells{cellnum}.smb_centers.go,fliplr(V1cells{cellnum}.smb_centers.go)],...
    'YData',[cellfun(@(x) prctile(x,5),randsmb80),fliplr(cellfun(@(x) prctile(x,95),randsmb80))],'FaceColor',[0 1 0],...
    'EdgeColor','none','FaceAlpha',0.1);
f=gca;
f.YTick = [-50 50];
