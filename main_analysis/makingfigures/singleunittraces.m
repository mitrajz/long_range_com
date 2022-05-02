%% 
% Do with pl20 an150 ms and smoothing 3
% - Errorbars are 95% bootstrapped confidence interval of the mean. 
% - Everything is only for go trials
% - In the influence-over-time plots, the error bars are when shuffling
% the id of ls vs bs trials and calculating % difference, and the main line
% is just smb

%% BEFORE RUNNING: You just need to load a cell
nrepbtstrp = 100; % def 1000
alphabtstrp = 5;
% which plots to do
smoothingspan = 3; % only used for traces (1: no smoothing)


% timeoffset determines if the point showa tha activity before, after, or
% middle of the bin: if timeoffset = 0; the shown point is acausal (shows
% the activity after it, if timeoffset = edgestep/(samplingf/1000); it is
% causal, and as default, it is halfway
timeoffset = .5*params.edgestep_pl/(params.samplingf/1000);

% plot appearance params
% yaxis hight ans xaxis
yl = 70;
xl = [-100 600];
% optopatch hight
OpatchH = 200;
% trace properties inside and outside of optopatch
maintrace.ls = '--';
maintrace.lw = 0.5;
maintrace.lc.bs = [0.2 0.2 0.2];
maintrace.lc.ls = [0 0.5 1];
maintrace.alpha = 0.1;

border = 10;

subtrace.ls = '-';
subtrace.lw = 1.5;
subtrace.lc.bs = [0.2 0.2 0.2];
subtrace.lc.ls = [0 0.5 1];
subtrace.alpha = 0.2;

%% raw traces og LM and V1 (all units) in go vs. nogo

xpoints = ((-params.Window):1*params.edgestep_pl:(params.Window-1))/(params.samplingf/1000);


targetcellname = 'LMcells{14}';
targetcell = eval(targetcellname);


%% traces
figure;
for lag = 1:8
    
    subplot(8,1,lag)
    
    % opto patch
    x = [targetcell.smb_centers.go(lag) targetcell.smb_centers.go(lag),...
        targetcell.smb_centers.go(lag)+150 targetcell.smb_centers.go(lag)+150];
    
    y = [0 OpatchH ...
        OpatchH 0];
    hold on;
    patch(x,y,[0 0 1],'FaceAlpha',0.1,'EdgeColor','none')
    %
    
    allbs = smooth(targetcell.nbsAv.go,smoothingspan)';
    allls = smooth(targetcell.nlsAv.go(lag,:),smoothingspan)';
    
    
    allbs_rand_mean=nan(nrepbtstrp,size(allbs,2));
    allls_rand_mean=nan(nrepbtstrp,size(allls,2));
    for randsamp = 1:nrepbtstrp
        allbs_rand =  smooth(nanmean(targetcell.nbs.go(randi(size(targetcell.nbs.go,1),1,size(targetcell.nbs.go,1)),:),1) , smoothingspan)';
        allls_rand = smooth(nanmean(targetcell.nls.go{lag}(randi(size(targetcell.nls.go{lag},1),1,size(targetcell.nls.go{lag},1)),:),1) , smoothingspan)';;
        allbs_rand_mean(randsamp,:) = nanmean(allbs_rand,1)/(params.edgestep_pl/30);
        allls_rand_mean(randsamp,:) = nanmean(allls_rand,1)/(params.edgestep_pl/30);
    end
    
    % main traces and btstrp CI
    %%% bs
    plot(timeoffset + xpoints,1000*nanmean(allbs,1)/(params.edgestep_pl/30),...
        'LineStyle',maintrace.ls,'Color',maintrace.lc.bs,'LineWidth',maintrace.lw);
    
    hold on;patch([timeoffset + xpoints,fliplr(timeoffset + xpoints)],...
        [1000*prctile(allbs_rand_mean,100 - (alphabtstrp/2)),fliplr(1000*prctile(allbs_rand_mean,alphabtstrp/2))],...
        maintrace.lc.bs,'FaceAlpha',maintrace.alpha,'EdgeAlpha',0);
    %%% ls
    hold on;plot(timeoffset + xpoints,1000*nanmean(allls,1)/(params.edgestep_pl/30),...
        'LineStyle',maintrace.ls,'Color',maintrace.lc.ls,'LineWidth',maintrace.lw);
    
    hold on;patch([timeoffset + xpoints,fliplr(timeoffset + xpoints)],...
        [1000*prctile(allls_rand_mean,100 - (alphabtstrp/2)),fliplr(1000*prctile(allls_rand_mean,alphabtstrp/2))],...
        maintrace.lc.bs,'FaceAlpha',maintrace.alpha,'EdgeAlpha',0);
    
    % patch traces
    allt = timeoffset + xpoints;
    patchind = find(allt>=x(1)-border & allt<=x(3)+border);
    % bs
    plot(allt(patchind),1000*nanmean(allbs(patchind),1)/(params.edgestep_pl/30),...
        'LineStyle',subtrace.ls,'Color',subtrace.lc.bs,'LineWidth',subtrace.lw);
    
    hold on;patch([allt(patchind),fliplr(allt(patchind))],...
        [1000*prctile(allbs_rand_mean(:,patchind),100 - (alphabtstrp/2)),fliplr(1000*prctile(allbs_rand_mean(:,patchind),alphabtstrp/2))],...
        subtrace.lc.bs,'FaceAlpha',subtrace.alpha,'EdgeAlpha',0);
    % ls
    hold on;plot(allt(patchind),1000*nanmean(allls(patchind),1)/(params.edgestep_pl/30),...
        'LineStyle',subtrace.ls,'Color',subtrace.lc.ls,'LineWidth',subtrace.lw);
    
    hold on;patch([allt(patchind),fliplr(allt(patchind))],...
        [1000*prctile(allls_rand_mean(:,patchind),100 - (alphabtstrp/2)),fliplr(1000*prctile(allls_rand_mean(:,patchind),alphabtstrp/2))],...
        subtrace.lc.ls,'FaceAlpha',subtrace.alpha,'EdgeAlpha',0);
    
    
    
    % figure appearance
    xlabel('Time (ms)');
    xlim(xl);ylim([0 yl]);
    
    axx = gca;
    % hold on; patch([0 stimlengthinms stimlengthinms 0],[axx.YLim(1) axx.YLim(1) axx.YLim(2) axx.YLim(2)],[0.7 0.7 0.7],...
    %     'FaceAlpha',0.3,'EdgeColor','none')
    axx.YTick=[];axx.XTick=[];
    if lag == 8
        axx.XTick=0:100:1000;
        axx.XTickLabel={'0','','200','','400','','600','','800','','1000'};
    end
    if lag ==1
        axx.YTick=0:yl/5:yl;
        axx.YTickLabel=arrayfun(@(x) num2str(x),axx.YTick,'UniformOutput',0);
        ylabel('Firing rate (Hz)');
    end
    set(gcf,'Color','w')
    box(axx,'off');
    
    set(axx,'Position',[0.13,(0.95-0.1*lag),0.775,0.1])
end

%%

figure;%plot(targetcell.smb_centers.go,targetcell.smb.go,'g.-')

difs = cellfun(@(x,y) 100*(x- nanmean(y))./nanmean(y),targetcell.laAls.go,targetcell.laAbs.go,'UniformOutput',0);
figure;errorbar(targetcell.smb_centers.go,cellfun(@(x) nanmean(x),difs),cellfun(@(x) 2*nanstd(x)/sqrt(numel(x)),difs),'g.-')

box('off')
ylabel('% Influence')
xlabel('Time(ms)')
ylim([-70 70])
xlim([-100 600])

%%% patch 1 with 95% ci of baseline: gray
hold on;patch('XData',[targetcell.smb_centers.go,fliplr(targetcell.smb_centers.go)],...
    'YData',[cellfun(@(x) 100*2*nanstd(x)./sqrt(length(x))./nanmean(x),targetcell.laAbs.go),...
    fliplr(cellfun(@(x) 100*-2*nanstd(x)./sqrt(length(x))./nanmean(x),targetcell.laAbs.go))],'FaceColor',[0 0 0],...
    'EdgeColor','none','FaceAlpha',0.1);


%%%% patch 2 with randomizing laser vs baseline: green
randsmb80 = cell(1,8);
for nlag=1:8
    randrep=1000;
    randsmb80{nlag} = nan(1,randrep);
    for i=1:randrep
        pool=[targetcell.laAls.go{nlag};targetcell.laAbs.go{nlag}];
        pool=pool(randperm(length(pool)));
        nlasertrials= size(targetcell.laAls.go{nlag},1);
        randsmb80{nlag}(i) =  100*(mean(pool(1:nlasertrials)) - mean(pool(nlasertrials+1:end)))/...
             mean(pool(nlasertrials+1:end));
        
    end
end


hold on;patch('XData',[targetcell.smb_centers.go,fliplr(targetcell.smb_centers.go)],...
    'YData',[cellfun(@(x) prctile(x,5),randsmb80),fliplr(cellfun(@(x) prctile(x,95),randsmb80))],'FaceColor',[0 1 0],...
    'EdgeColor','none','FaceAlpha',0.1);
f=gca;
f.YTick = [-50 -25 0  25 50];
