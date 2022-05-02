%% IMPORTANT: this code is an externsion of singleunittraces.m to allunits instead of 1
% For information read documentation of singleunittraces.m

nrepbtstrp = 10; % def 100, set to 10 for faster run
alphabtstrp = 5;
normalize = 0;
% which plots to do
smoothingspan = 1; % only used for traces


% timeoffset determines if the point showa tha activity before, after, or
% middle of the bin: if timeoffset = 0; the shown point is acausal (shows
% the activity after it, if timeoffset = edgestep/(samplingf/1000); it is
% causal, and as default, it is halfway
timeoffset = .5*params.edgestep_pl/(params.samplingf/1000);

% plot appearance params
% yaxis hight ans xaxis
yl = 18; % 18
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

%subselect_ind_phrase = sprintf('find(cellfun(@(x) x.spikewidth,%s)>th & cellfun(@(x) x.spikewidth,%s)<th2)','V1cells','V1cells')
%subselect_ind_phrase = sprintf('find(cellfun(@(x) x.spikewidth,%s)<=th)','V1cells');
%targetcellname = 'LMcells(eval(subselect_ind_phrase))';
targetcellname = 'V1cells';

targetcell = eval(targetcellname);


%% traces
figure;
for lag = 1:8
    
    subplot(8,1,lag)
    
    % opto patch
    x = [targetcell{1}.smb_centers.go(lag) targetcell{1}.smb_centers.go(lag),...
        targetcell{1}.smb_centers.go(lag)+150 targetcell{1}.smb_centers.go(lag)+150];
    
    y = [0 OpatchH ...
        OpatchH 0];
    hold on;
    patch(x,y,[0 0 1],'FaceAlpha',0.1,'EdgeColor','none')
    %
    
    
    allbs=cell2mat(transpose(cellfun(@(x) smooth(x.nbsAv.go,smoothingspan)' ,targetcell,'UniformOutput',0)));
    allls=cell2mat(transpose(cellfun(@(x) smooth(x.nlsAv.go(lag,:),smoothingspan)' ,targetcell,'UniformOutput',0)));
    
    if normalize
        allbs = allbs./repmat(cellfun(@(x) mean(x.stimspikes.go),targetcell)',1,size(allbs,2));
        allls = allls./repmat(cellfun(@(x) mean(x.stimspikes.go),targetcell)',1,size(allls,2));
    end
    
    allbs_rand_mean=nan(nrepbtstrp,size(allbs,2));
    allls_rand_mean=nan(nrepbtstrp,size(allls,2));
    for randsamp = 1:nrepbtstrp
        allbs_rand = cell2mat(transpose(cellfun(@(x) smooth(nanmean(x.nbs.go(randi(size(x.nbs.go,1),1,size(x.nbs.go,1)),:),1) , smoothingspan)',...
            targetcell,'UniformOutput',0)));
        allls_rand = cell2mat(transpose(cellfun(@(x) smooth(nanmean(x.nls.go{lag}(randi(size(x.nls.go{lag},1),1,size(x.nls.go{lag},1)),:),1) , smoothingspan)',...
            targetcell,'UniformOutput',0)));
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
    plot(allt(patchind),1000*nanmean(allbs(:,patchind),1)/(params.edgestep_pl/30),...
        'LineStyle',subtrace.ls,'Color',subtrace.lc.bs,'LineWidth',subtrace.lw);
    
    hold on;patch([allt(patchind),fliplr(allt(:,patchind))],...
        [1000*prctile(allbs_rand_mean(:,patchind),100 - (alphabtstrp/2)),fliplr(1000*prctile(allbs_rand_mean(:,patchind),alphabtstrp/2))],...
        subtrace.lc.bs,'FaceAlpha',subtrace.alpha,'EdgeAlpha',0);
    % ls
    hold on;plot(allt(patchind),1000*nanmean(allls(:,patchind),1)/(params.edgestep_pl/30),...
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



