% summary of the approaches:
% 1 - combined pvalue and effect size: if the min pvalue and max effect size
% happen at the same lag, that lag is considered effect onset. Advantage:
% it is laser aligned and fine precision (lagstepsize) inspite of the
% larger window of spike counting for effect(lagwinsize). Disadvantage:
% doesn't detect early effects, but more like max effect, that can happen
% much later. Especially if baseline activity increases, while the laser
% activity is down to zero (due to this problem, it is probably more
% reasonable to average over lags for each cell)
%
% [minpval , minpvalindex]=min(smooth(targetcell{i}.pval_lags.go{l},smoothparam));
% [maxeffect , maxeffectindex]=max(smooth(targetcell{i}.effectsize_lags.go{l},smoothparam));
% if minpvalindex == maxeffectindex
%    onset_lag_ms.go{i}(l) = (minpvalindex-1)*params.movingwin_withlags.lagsize/30;
% end
%
% 2 - purely based on pvalue: The first time pvalue passes 0.05 threshold.
% Seems to work fine for strong effects (initial site of inhibition), but
% for the other area, not so well
%
% minpvalindex = find(smooth(targetcell{i}.pval_lags.go{l},smoothparam)<0.05,1);
% if numel(minpvalindex)
%   onset_lag_ms.go{i}(l) = (minpvalindex-1)*params.movingwin_withlags.lagsize/30;
% end
%
% 3 - Purely based on effect size: The first time effect size crosses a
% certain threshold. The threshold is based on zscore of effectsize pre
% laser (might introduce biad for different lags). The advantage over
% pvalue is that the variability is based on laser vs. baseline differences
% pre laser, instead of the separability of laser vs. baseline trials after
% So, we get the time post laser that laser and baseline are significantly
% different from each other, not with respect to the distribution of each
% condition, but with respect to how different they were pre laser.
% -- In general lagged window with respect to laser onset is preferred to
% binned averages. Due to (1) laser alignement, (2) fine resolution of
% lags, independent of the size of window, which is chosen much larger to
% integrate more spikes


%%
% load(cellfiletoload);
%
% [LMcells,V1cells,Behcells,params] = cleanupcells(cleanupcriteria,LMcells,V1cells,Behcells,params);
% [cell_keep_ind,LMcells,V1cells,Behcells,params] = removecells(cell_keep_ind,cellremovecriteria,LMcells,V1cells,Behcells,params);

%% params

% on/offset calculation parameters
smoothparam = 1;
zscorethreshold  = 2;
onsetlimit = 150*30/params.movingwin_withlags.lagsize; % 150ms: duration of laser
% histogram plots
hist_yl = [0 0.02]; % ylimit
hist_edges = 0:6:params.movingwin_withlags.nlags * params.movingwin_withlags.lagsize/30; % covered time in ms
% plotting of example cells
nexamples = 4;
samplingf = 30000;
Window = 30000;
stecoeff = 1;
inhpatchdur_ms =150;
edgestep_pl = params.edgestep_pl;
timeoffset = .5*edgestep_pl/(samplingf/1000);
onsetplotsteps = 4;

%%

% get onset and offset based on the first and last point to have high
% zscores compared to pre laser baseline

targetcell = eval(targetcellname);
% subselect cells for plots only
if subselectcells
    subselect_ind = eval(subselect_ind_phrase);
    fprintf('subselecting %d cells\n',numel(subselect_ind))
    targetcell = targetcell(subselect_ind);
    if showspikeshape
        figure;plot(allspikes(:,subselect_ind),'k')
    end
end

onset_lag_ms.go = cell(1,length(targetcell));
onset_lag_ms.nogo = cell(1,length(targetcell));
offset_lag_ms.go = cell(1,length(targetcell));
offset_lag_ms.nogo = cell(1,length(targetcell));

for i = 1:length(targetcell)
    % go
    onset_lag_ms.go{i} = nan(1,length(targetcell{i}.effectsize_lags.go));
    offset_lag_ms.go{i} = nan(1,length(targetcell{i}.effectsize_lags.go));
    for l =1:length(targetcell{i}.effectsize_lags.go)
        
        zscoredeffect = (targetcell{i}.effectsize_lags.go{l}-nanmean(targetcell{i}.effectsize_lags_pre.go{l}))/nanstd(targetcell{i}.effectsize_lags_pre.go{l});
        minpvalindex_first = find(abs(zscoredeffect)>=zscorethreshold,1,'first');
        minpvalindex_last = find(abs(zscoredeffect)>=zscorethreshold,1,'last');
        % if a point passing the zscore threshold is detected, is inside
        % the allowed zone,and at least one point later, it is still above
        % threshold
        if numel(minpvalindex_first) & minpvalindex_first<min(onsetlimit,length(zscoredeffect)) & abs(zscoredeffect(minpvalindex_first+1))>=zscorethreshold
            onset_lag_ms.go{i}(l) = (minpvalindex_first-1)*params.movingwin_withlags.lagsize/30;
            if minpvalindex_first ~= minpvalindex_last
                offset_lag_ms.go{i}(l) = (minpvalindex_last-1)*params.movingwin_withlags.lagsize/30;
            end
        end
        
    end
    % nogo
    onset_lag_ms.nogo{i} = nan(1,length(targetcell{i}.effectsize_lags.nogo));
    offset_lag_ms.nogo{i} = nan(1,length(targetcell{i}.effectsize_lags.nogo));
    for l = 1:length(targetcell{i}.effectsize_lags.nogo)
        
        zscoredeffect = (targetcell{i}.effectsize_lags.nogo{l}-nanmean(targetcell{i}.effectsize_lags_pre.nogo{l}))/nanstd(targetcell{i}.effectsize_lags_pre.nogo{l});
        minpvalindex_first = find(abs(zscoredeffect)>=zscorethreshold,1,'first');
        minpvalindex_last = find(abs(zscoredeffect)>=zscorethreshold,1,'last');
        
        if numel(minpvalindex_first) & minpvalindex_first<min(onsetlimit,length(zscoredeffect)) & abs(zscoredeffect(minpvalindex_first+1))>=zscorethreshold
            onset_lag_ms.nogo{i}(l) = (minpvalindex_first-1)*params.movingwin_withlags.lagsize/30;
            if minpvalindex_first ~= minpvalindex_last
                offset_lag_ms.nogo{i}(l) = (minpvalindex_last-1)*params.movingwin_withlags.lagsize/30;
            end
        end
    end
    
end
effect_lag_ms = eval([effect,'_lag_ms']);

%% 2 types of plots with go/nogo separation
% histogram of each cell, averaged over lags
figure;
s1 = subplot(2,1,1);
gomed = nanmedian(nanmean(cell2mat(effect_lag_ms.go'),2));
nogomed = nanmedian(nanmean(cell2mat(effect_lag_ms.nogo'),2));
hold on;histogram(nanmean(cell2mat(effect_lag_ms.go'),2),hist_edges,'Normalization',globalplot.histnorm,'FaceColor','g','FaceAlpha',0.1,'EdgeAlpha',0);
hold on;histogram(nanmean(cell2mat(effect_lag_ms.go'),2),hist_edges,'Normalization',globalplot.histnorm,'EdgeColor',[0 0.8 0],'EdgeAlpha',1,'DisplayStyle','stairs');
hold on;line([gomed,gomed],hist_yl,'Color',[0 0.9 0])
hold on;histogram(nanmean(cell2mat(effect_lag_ms.nogo'),2),hist_edges,'Normalization',globalplot.histnorm,'FaceColor','r','FaceAlpha',0.1,'EdgeAlpha',0);
hold on;histogram(nanmean(cell2mat(effect_lag_ms.nogo'),2),hist_edges,'Normalization',globalplot.histnorm,'EdgeColor',[0.8 0 0],'EdgeAlpha',1,'DisplayStyle','stairs');
hold on;line([nogomed,nogomed],hist_yl,'Color',[0.9 0 0])
text([max(nogomed,gomed)],.9*hist_yl(2),sprintf('p = %f',ranksum(nanmean(cell2mat(effect_lag_ms.go'),2),nanmean(cell2mat(effect_lag_ms.nogo'),2))),'FontSize',6)
s1.Title.String = sprintf('onset of change for cells in %s, averaged over all silencing onsets for each cell',targetcellname(1:2));
s1.Title.FontSize = 8;
s1.Title.FontWeight = 'normal';
% histogram of each cell all lags included
s2 = subplot(2,1,2);
gomed = nanmedian(reshape(cell2mat(effect_lag_ms.go'),1,[]));
nogomed = nanmedian(reshape(cell2mat(effect_lag_ms.nogo'),1,[]));
hold on; histogram(reshape(cell2mat(effect_lag_ms.go'),1,[]),hist_edges,'Normalization',globalplot.histnorm,'FaceColor','g','FaceAlpha',0.1,'EdgeAlpha',0);
hold on; histogram(reshape(cell2mat(effect_lag_ms.go'),1,[]),hist_edges,'Normalization',globalplot.histnorm,'EdgeColor',[0 0.8 0],'EdgeAlpha',1,'DisplayStyle','stairs');
hold on;line([gomed,gomed],hist_yl,'Color',[0 0.9 0])
hold on;histogram(reshape(cell2mat(effect_lag_ms.nogo'),1,[]),hist_edges,'Normalization',globalplot.histnorm,'FaceColor','r','FaceAlpha',0.1,'EdgeAlpha',0);
hold on;histogram(reshape(cell2mat(effect_lag_ms.nogo'),1,[]),hist_edges,'Normalization',globalplot.histnorm,'EdgeColor',[0.8 0 0],'EdgeAlpha',1,'DisplayStyle','stairs');
hold on;line([nogomed,nogomed],hist_yl,'Color',[0.9 0 0])
text([max(nogomed,gomed)],.9*hist_yl(2),sprintf('p = %f',ranksum(reshape(cell2mat(effect_lag_ms.go'),1,[]),reshape(cell2mat(effect_lag_ms.nogo'),1,[]))),'FontSize',6)
s2.Title.String = sprintf('onset of change for cells in %s, all silencing onsets and cells combined',targetcellname(1:2));
s2.Title.FontSize = 8;
s2.Title.FontWeight = 'normal';
% separate for different lags
figure;
allgoeffects = cell2mat(effect_lag_ms.go');
allnogoeffects = cell2mat(effect_lag_ms.nogo');
for l=1:length(targetcell{i}.effectsize_lags.go)
    
    s = subplot(length(targetcell{i}.effectsize_lags.go),1,l);
    gomed = nanmedian(nanmean(allgoeffects(:,l),2));
    nogomed = nanmedian(nanmean(allnogoeffects(:,l),2));
    hold on;histogram(nanmean(allgoeffects(:,l),2),hist_edges,'Normalization',globalplot.histnorm,'FaceColor','g','FaceAlpha',0.1,'EdgeAlpha',0);
    hold on;histogram(nanmean(allgoeffects(:,l),2),hist_edges,'Normalization',globalplot.histnorm,'EdgeColor',[0 0.8 0],'EdgeAlpha',1,'DisplayStyle','stairs');
    hold on;line([gomed,gomed],hist_yl,'Color',[0 0.9 0])
    hold on;histogram(nanmean(allnogoeffects(:,l),2),hist_edges,'Normalization',globalplot.histnorm,'FaceColor','r','FaceAlpha',0.1,'EdgeAlpha',0);
    hold on;histogram(nanmean(allnogoeffects(:,l),2),hist_edges,'Normalization',globalplot.histnorm,'EdgeColor',[0.8 0 0],'EdgeAlpha',1,'DisplayStyle','stairs');
    hold on;line([nogomed,nogomed],hist_yl,'Color',[0.9 0 0])
    text([max(nogomed,gomed)],.9*hist_yl(2),sprintf('p = %f',ranksum(nanmean(allgoeffects(:,l),2),nanmean(allnogoeffects(:,l),2))),'FontSize',6)
    s=gca;
    if l==1
        s.Title.String = sprintf('onset of change for cells in %s, averaged over all silencing onsets for each cell',targetcellname(1:2));
        s.Title.FontSize = 8;
        s.Title.FontWeight = 'normal';
    end
    xlabel(sprintf('silencing onset %d',l))
end
% separate for different lags
%% 2 types of plots, combining go and nogo
% go and nogo are combined as different cells
effect_lag_ms_both = [cell2mat(effect_lag_ms.go');cell2mat(effect_lag_ms.nogo')];

figure;
s1 = subplot(2,1,1);
med = nanmedian(nanmean(effect_lag_ms_both,2));
hold on;histogram(nanmean(effect_lag_ms_both,2),hist_edges,'Normalization',globalplot.histnorm,'FaceColor','k','FaceAlpha',0.9,'EdgeAlpha',0);
hold on;line([med,med],hist_yl,'Color',[0 0 0])
text([med+1],.9*hist_yl(2),sprintf('%d ms',round(med)))
s1.Title.String = sprintf('onset of change for cells in %s, averaged over all silencing onsets for each cell',targetcellname(1:2));
s1.Title.FontSize = 8;
s1.Title.FontWeight = 'normal';

s2 = subplot(2,1,2);
med = nanmedian(reshape(effect_lag_ms_both,1,[]));
hold on; histogram(reshape(effect_lag_ms_both,1,[]),hist_edges,'Normalization',globalplot.histnorm,'FaceColor','k','FaceAlpha',0.9,'EdgeAlpha',0);
hold on;line([med,med],hist_yl,'Color',[0 0 0])
text([med+1],.9*hist_yl(2),sprintf('%d ms',round(med)))
s2.Title.String = sprintf('onset of change for cells in %s, all silencing onsets and cells combined',targetcellname(1:2));
s2.Title.FontSize = 8;
s2.Title.FontWeight = 'normal';
% separate for different lags
figure;
for l=1:length(targetcell{i}.effectsize_lags.go)
    
    s1 = subplot(length(targetcell{i}.effectsize_lags.go),1,l);
    med = nanmedian(nanmean(effect_lag_ms_both(:,l),2));
    hold on;histogram(nanmean(effect_lag_ms_both(:,l),2),hist_edges,'Normalization',globalplot.histnorm,'FaceColor','k','FaceAlpha',0.9,'EdgeAlpha',0);
    hold on;line([med,med],hist_yl,'Color',[0 0 0])
    text([med+1],.9*hist_yl(2),sprintf('%d ms',round(med)))
    
    if l==1
        s.Title.String = sprintf('onset of change for cells in %s, averaged over all silencing onsets for each cell',targetcellname(1:2));
        s.Title.FontSize = 8;
        s.Title.FontWeight = 'normal';
    end
    xlabel(sprintf('silencing onset %d',l))
end
%% Example traces, here only for go conditions to check the performance of onset/offset detection
figure;
cellagtable = cell2mat(effect_lag_ms.go');
allags = unique(cellagtable(~isnan(cellagtable)));
for i = 1:onsetplotsteps:length(allags)
    % for each value of onset/offset find cellid and lags corresponding to
    % id=t. row:cell, cul:lag - The same number of rows and culs - and they
    % should be jointly manipulated or indexed
    [row,cul] = find(cellagtable == allags(i));
    randorder = randperm(size(row,1));
    % at most nexamples plots per row
    for i_cellid = 1:min(nexamples,size(row,1))
        % choose the subplot
        s = subplot(ceil(length(allags)/onsetplotsteps),nexamples,(i-1)/onsetplotsteps*nexamples + i_cellid);
        % choose the example (cell. lag)
        % the same index used for both row and cul
        cellid = row(randorder(i_cellid));
        lag = cul(randorder(i_cellid));
        
        if false
            % another way to choose: with replacement
            randentryind = randi(length(cul));
            cellid = row(randentryind);
            lag = cul(randentryind);
        end
        
        % make basline and laser plots
        p=shadedErrorBar(timeoffset+((-Window):1*edgestep_pl:(Window-1))/(samplingf/1000),...
            1000*nanmean(targetcell{cellid}.nbs.go,1)/(edgestep_pl/30),...
            1000*stecoeff*nanstd(targetcell{cellid}.nbs.go,1)/(edgestep_pl/30)/sqrt(size(targetcell{cellid}.nbs.go,1)),...
            {'Color',[0 0 0],'LineWidth',1.6},1);
        hold on;
        shadedErrorBar(timeoffset+((-Window):1*edgestep_pl:(Window-1))/(samplingf/1000),...
            1000*nanmean(targetcell{cellid}.nls.go{lag},1)/(edgestep_pl/30),...
            1000*stecoeff*nanstd(targetcell{cellid}.nls.go{lag},1)/(edgestep_pl/30)/sqrt(size(targetcell{cellid}.nls.go{lag},1)),...
            {'Color',[0 0 1],'LineWidth',1.6},1);
        
        % mark laser duration, onset and offset points, labels, set axis
        % limit
        startpoint=targetcell{cellid}.smb_centers.go(lag);
        yl = get(gca,'YLim');
        hold on;patch([ startpoint startpoint+inhpatchdur_ms startpoint+inhpatchdur_ms startpoint],[0 0 yl(2) yl(2)],'b','FaceAlpha',0.1,'EdgeAlpha',0)
        hold on;line([startpoint+onset_lag_ms.go{cellid}(lag),startpoint+onset_lag_ms.go{cellid}(lag)],yl,'LineStyle','--','Color','k')
        hold on;line([startpoint+offset_lag_ms.go{cellid}(lag),startpoint+offset_lag_ms.go{cellid}(lag)],yl,'LineStyle','--','Color','r')
        xlim([startpoint-200 startpoint+inhpatchdur_ms+300])
        if i_cellid == 1
            s.YLabel.String = sprintf('%s =%d ms',effect,allags(i));
            s.YLabel.FontSize = 7;
        end
        % xlabel('Time (ms)');ylabel('Firing rate (Hz)');
    end
    
end
sprintf('found %s values for %f %s of cases in %s (rest nan)',effect,100*numel(find(~isnan(cellagtable)))/numel(find(cellagtable)),'%',targetcellname(1:2))


