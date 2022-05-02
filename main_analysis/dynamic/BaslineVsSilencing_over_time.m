% This scripts plots activity with/out laser either in all cells or in single cells in a given are
% and when used for single cells, also plots the effcet over time and shuffle chance level

%% load file, remove and clear cells (first part of aggregate analysis)

%% params: set/check all before running

stimlengthinms = 500; 
nrepbtstrp = 100;
alphabtstrp = 5;
mean_standarderr_coeff = 1.96; % for the effect plots

% select what to plot
areaind = 1; % 1:V1, 2:LM
wide_spiking = 0;
selectcell =53; % 
singlecellylim = 100;
smoothingspan = 3;%3
trialtype = 'nogo';



lag=2;
alllags = 1;
avgonogo = 0; % only to be used with untrained data, makes extra virtual cells with nogo and nogo fields switched

%%

timeoffset = .5*params.edgestep_pl/(params.samplingf/1000);
plotoffset = 0;
lagrange=lag;
if alllags
    clear lag
    lagrange= 1:8;
end

xpoints = ((-params.Window):1*params.edgestep_pl:(params.Window-1))/(params.samplingf/1000);
targetcellnamelist = {'V1cells','LMcells'};


laserXpoints = nanmean(unique( [cell2mat(cellfun(@(x) x.smb_centers.go,V1cells,'UniformOutput',0)');...
    cell2mat(cellfun(@(x) x.smb_centers.go,LMcells,'UniformOutput',0)')],...
    'rows'),1);


%% making targetcell
targetcellname = targetcellnamelist{areaind};
targetcell = eval(targetcellname);
if avgonogo
    targetcell = virtualgonocomb(targetcell);
end

if  selectcell > 0
    targetcell = targetcell( selectcell);
end
if wide_spiking
    targetcell = targetcell(find(cellfun(@(x) x.spikewidth,targetcell)>th & cellfun(@(x) x.spikewidth,targetcell)<th2));
end

%% choosing the right trial type
for n=1:length(targetcell)
    fn=fieldnames(targetcell{n});
    for f=1:length(fn)
        if isfield(getfield(targetcell{n},fn{f}),trialtype)
            targetfield = getfield(getfield(targetcell{n},fn{f}),trialtype);
            targetcell{n} = setfield(targetcell{n},fn{f},targetfield);
        end
    end
end
%%

figure;

for lag = lagrange
    
    allbs=cell2mat(transpose(cellfun(@(x) smooth(x.nbsAv,smoothingspan)' ,targetcell,'UniformOutput',0)));
    allls=cell2mat(transpose(cellfun(@(x) smooth(x.nlsAv(lag,:),smoothingspan)' ,targetcell,'UniformOutput',0)));
    
    %%% traces
    if plotstodo.trace
        plot(timeoffset + xpoints,plotoffset+1000*nanmean(allbs,1)/(params.edgestep_pl/30),'Color',[0 0 0],'LineWidth',1.6);
        hold on;plot(timeoffset + xpoints,plotoffset+1000*nanmean(allls,1)/(params.edgestep_pl/30),'Color',[0 0 1],'LineWidth',1.6);
        %%% errorbars by bootstrapping across trials
        if doerbars
            allbs_rand_mean=nan(nrepbtstrp,size(allbs,2));
            allls_rand_mean=nan(nrepbtstrp,size(allls,2));
            for randsamp = 1:nrepbtstrp
                allbs_rand = cell2mat(transpose(cellfun(@(x) smooth(nanmean(x.nbs(randi(size(x.nbs,1),1,size(x.nbs,1)),:)) , smoothingspan)',...
                    targetcell,'UniformOutput',0)));
                allls_rand = cell2mat(transpose(cellfun(@(x) smooth(nanmean(x.nls{lag}(randi(size(x.nls{lag},1),1,size(x.nls{lag},1)),:)) , smoothingspan)',...
                    targetcell,'UniformOutput',0)));
                allbs_rand_mean(randsamp,:) = nanmean(allbs_rand,1)/(params.edgestep_pl/30);
                allls_rand_mean(randsamp,:) = nanmean(allls_rand,1)/(params.edgestep_pl/30);
            end
            hold on;patch([timeoffset + xpoints,fliplr(timeoffset + xpoints)],...
                plotoffset+[1000*prctile(allbs_rand_mean,100 - (alphabtstrp/2)),fliplr(1000*prctile(allbs_rand_mean,(alphabtstrp/2)))],[0 0 0],'FaceAlpha',0.1,'EdgeAlpha',0);
            hold on;patch([timeoffset + xpoints,fliplr(timeoffset + xpoints)],...
                plotoffset+[1000*prctile(allls_rand_mean,100 - (alphabtstrp/2)),fliplr(1000*prctile(allls_rand_mean,(alphabtstrp/2)))],[0 0 1],'FaceAlpha',0.1,'EdgeAlpha',0);
                        
        end
        
        %%% plot properties and patches
        
        xlabel('Time (ms)');ylabel('Firing rate (Hz)');
        xlim([-100 1000]);
        ylim([0 25]);
        if selectcell
            ylim([0  singlecellylim])
        end        
        axx = gca;        
        startpoint=targetcell{1}.smb_centers(lag);
        dur3=150;
        hold on;patch([ startpoint startpoint+dur3 startpoint+dur3 startpoint],...
            [plotoffset+0 plotoffset+0 plotoffset+axx.YLim(2) plotoffset+axx.YLim(2)],'b','FaceAlpha',0.1,'EdgeAlpha',0)
        
        hold on; patch([0 stimlengthinms stimlengthinms 0],[axx.YLim(2)*.97 axx.YLim(2)*.97 axx.YLim(2) axx.YLim(2)],[0.1 0.1 0.1],...
            'FaceAlpha',0.3,'EdgeColor','none')
        axx.XTick=0:100:1000;
        axx.XTickLabel={'0','','200','','400','','600','','800','','1000'};
        
        axx.YTick=0:5:20;
        axx.YTickLabel={'0','5','10','15','20','25'};
        if selectcell
            axx.YTick=0:10:singlecellylim;
            axx.YTickLabel=arrayfun(@(x) num2str(x), 0:10:singlecellylim,'UniformOutput',0);
        end
        set(gcf,'Color','w')
        legend({sprintf('%s - correct - baseline',trialtype),sprintf('%s - correct - silenced',trialtype)})
        legend('boxoff')
        box(axx,'off');
        %%% 
    end
    
    if alllags
        plotoffset  =  plotoffset - axx.YLim(2);
    end
end
if alllags
    axx.YLim  =  [plotoffset+axx.YLim(2)  axx.YLim(2)];
end
xlim([-100 600])

%% effect ovetime for single cells + chance level
colorcode = [0 0 0];
if strcmp(trialtype,'go')
    colorcode = [0 0.8 0];
elseif strcmp(trialtype,'nogo')
    colorcode = [0.8 0 0];
end

if selectcell>0 && alllags
    i=1 % because cell has been already selected
    chancelevel = nan(1,8);
    for l=1:8
        mean_standarderr = nanstd(targetcell{i}.laAbs{l})/sqrt(numel(targetcell{i}.laAbs{l}));
        chancelevel(l) = 100 * mean_standarderr_coeff * mean_standarderr/targetcell{i}.smb_bs(l);
    end
    figure;plot(laserXpoints,targetcell{i}.smb,...
        'Color',colorcode,'LineWidth',1.5,'Marker','.','MarkerSize',10);
    hold on;
    patch([laserXpoints,fliplr(laserXpoints)],...
        [-chancelevel,fliplr(chancelevel)],colorcode,'EdgeAlpha',0,'FaceAlpha',0.1)
end