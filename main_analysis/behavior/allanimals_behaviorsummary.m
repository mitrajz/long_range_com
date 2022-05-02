


animallist = {'MPV17','MPV18_2',...
    'VL53','VL52','VL51','VL66','MPV35_2',...
    'VL61','VL63','VL55','VL59','VL50',...
    'MPV32_2','MPV33','MPV31','MPV34_2'};
preprocessinglist = {'2020_03_02_13_10_37','2020_03_02_13_27_45',...
    '2017_12_05_18_58_17','2017_12_04_15_08_47','2017_12_03_17_30_15','2018_03_20_11_53_19','2020_03_02_14_31_53',...
    '2018_03_16_16_53_13','2018_04_01_12_24_38','2018_04_01_15_03_48','2018_04_01_13_24_48','2017_12_01_16_21_49',...
    '2020_03_02_19_50_11','2020_03_02_19_35_12','2020_03_02_19_58_50','2020_03_02_20_32_35'};
exptype = {'FB','FB',...
    'FB','FB','FB','FB','FB',...
    'FF','FF','FF','FF','FF',...
    'FF','FF','FF','FF'};


dp = nan(1,length(animallist));
dp_nogrooming = nan(1,length(animallist));
H_nogrooming = nan(1,length(animallist));
F_nogrooming = nan(1,length(animallist));


allhitrate = cell(1,length(animallist));
allmissrate = cell(1,length(animallist));
allxpoints = cell(1,length(animallist));
allcrrate = cell(1,length(animallist));
allfarate = cell(1,length(animallist));
lickcumsum = cell(1,length(animallist));
licklatencydist = cell(1,length(animallist));

for animal_i=1:length(animallist)
    cd('/mnt/data/Mitra/figs');
    cd(sprintf('%s/preprocessing/%s',animallist{animal_i},preprocessinglist{animal_i}))
    matfilename = dir('task_*.mat');
    load(matfilename.name)
    perfmatfilename = dir('rev_perf*.mat');
    if ~isempty(perfmatfilename)
        load(perfmatfilename.name)
        % get first lick sample
        firstlicksample = nan(1,size(PAllOn,1));
        for trialind=1:1:size(PAllOn,1)
            fl=find(single(licks(trialind,floor(size(licks,2)/2):end)) == 1,1);
            if numel(fl) && fl>0
                firstlicksample(trialind) = fl;
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%% d-prime
    % this doesn't count early lick trials
    H = length(correctgotrialind)/(length(correctgotrialind)+length(incorrectgotrialind));
    
    F = length(incorrectnogotrialind)/(length(correctnogotrialind)+length(incorrectnogotrialind));
    
    dp(animal_i) = norminv(H) - norminv(F);
    
    % this doesn't count early lick trials
    H_nogrooming(animal_i) = length(intersect(correctgotrialind,nogroomingind))/...
        (length(intersect(correctgotrialind,nogroomingind))+length(intersect(incorrectgotrialind,nogroomingind)));
    
    F_nogrooming(animal_i) = length(intersect(incorrectnogotrialind,nogroomingind))/...
        (length(intersect(correctnogotrialind,nogroomingind))+length(intersect(incorrectnogotrialind,nogroomingind)));
    
    dp_nogrooming(animal_i) = norminv(H_nogrooming(animal_i)) - norminv(F_nogrooming(animal_i));
    
    %%%%%%%%%%%%%%%%%%%%% response latency and hit/fa/miss/cr rates, during both baseline and laser trials
    
    min_n_trials_per_delay = 25; 
    cumsumstep = 50; 
    onlynogrooming = 1;
    onlystableperiod = 0;
    normalizeforcorrect = 0;
    put2nanoutsidewin = 0;
    respwindowms = [100 650]; 
    onlycorrect = 0;
    
    
    [allhitrate{animal_i},allmissrate{animal_i},allfarate{animal_i},allcrrate{animal_i},allxpoints{animal_i}] = ...
        behsummaryhelper_rates(gotrialind,nogotrialind,...
        correctgotrialind,incorrectgotrialind,correctnogotrialind,incorrectnogotrialind,...
        stabletrialind,nogroomingind,LaserDelayBinned,...
        min_n_trials_per_delay,onlynogrooming,onlystableperiod,centers);
    
    
    [lickcumsum{animal_i},licklatencydist{animal_i}] = behsummaryhelper_lickhist(gotrialind,nogotrialind,...
        correctgotrialind,incorrectgotrialind,correctnogotrialind,incorrectnogotrialind,...
        stabletrialind,nogroomingind,LaserDelayBinned,...
        min_n_trials_per_delay,onlynogrooming,onlystableperiod,...
        cumsumstep,normalizeforcorrect,respwindowms,put2nanoutsidewin,centers,firstlicksample,onlycorrect);
    
    [chanceperf{animal_i}] = behsummaryhelper_chanceperf(gotrialind,nogotrialind,...
        correctgotrialind,incorrectgotrialind,correctnogotrialind,incorrectnogotrialind,...
        nogroomingind,onlynogrooming);
    
end
%% check xpoints
figure;plot(cell2mat(allxpoints')')
%% assess better than chance behavior: 99 percentile
figure;
plot_bar_width = 0.3;
for i=1:length(animallist)
    hold on; patch([i-plot_bar_width i+plot_bar_width i+plot_bar_width i-plot_bar_width],...
        [prctile(chanceperf{i},99.5) prctile(chanceperf{i},99.5),prctile(chanceperf{i},.5),prctile(chanceperf{i},.5)],...
        'k','FaceAlpha',0.1,'EdgeAlpha',0)
    hold on;scatter(i,dp_nogrooming(i),100,'k.');
end
set(gca,'XTick',1:length(animallist))
set(gca,'XTickLabel',animallist)
set(gca,'XTickLabelRotation',90)
ylabel('d-prime')

ind = find(dp_nogrooming>=cellfun(@(x) prctile(x,99.5),chanceperf));
%% plots
%%%%%%% params
ci_coef = 1.96;
%% basline plots - both FF and FB
%%%%%%%% baseline d prime excluding animals worse than thresold: for both FF and FB


figure;plot(dp_nogrooming,'.')

figure;ax0=gca;
boxplot(ax0,dp_nogrooming(ind),'Position',zeros(size(dp_nogrooming(ind))),...
    'PlotStyle','compact','MedianStyle','line','Colors',[0.7 0.7 0.7])
hold on;scatter(ax0,(rand(size(dp_nogrooming(ind)))*0.05)-0.025,dp_nogrooming(ind),100,'k.');ylim([0 3]);xlim([-1,1])
ylabel('d-prime')

figure;ax1=gca;
boxplot(ax1,H_nogrooming(ind),'Position',zeros(size(H_nogrooming(ind))),...
    'PlotStyle','compact','MedianStyle','line','Colors',[0.7 0.7 0.7])
hold on;scatter(ax1,(rand(size(H_nogrooming(ind)))*0.05)-0.025,H_nogrooming(ind),100,'k.');ylim([0 1]);xlim([-1,1])
ylabel('hit rate')

figure;ax2=gca;
boxplot(ax2,F_nogrooming(ind),'Position',zeros(size(F_nogrooming(ind))),...
    'PlotStyle','compact','MedianStyle','line','Colors',[0.7 0.7 0.7])
hold on;scatter(ax2,(rand(size(F_nogrooming(ind)))*0.05)-0.025,F_nogrooming(ind),100,'k.');ylim([0 1]);xlim([-1,1])
ylabel('false alarm rate')
%%%%%%% separate for FF and FB:
FFind = find(strcmp(exptype,'FF'));
FBind = find(strcmp(exptype,'FB'));
figure;ax=gca;
boxplot(ax,dp_nogrooming(intersect(ind,FFind)),'Position',zeros(size(dp_nogrooming(intersect(ind,FFind)))),...
    'PlotStyle','compact','MedianStyle','line','Colors',[0.7 0.7 0.7])
hold on;scatter(ax,(rand(size(dp_nogrooming(intersect(ind,FFind))))*0.05)-0.025,dp_nogrooming(intersect(ind,FFind)),100,'k.');
hold on;
boxplot(ax,dp_nogrooming(intersect(ind,FBind)),'Position',ones(size(dp_nogrooming(intersect(ind,FBind)))),...
    'PlotStyle','compact','MedianStyle','line','Colors',[0.7 0.7 0.7])
hold on;scatter(ax,(rand(size(dp_nogrooming(intersect(ind,FBind))))*0.05)-0.025+1,dp_nogrooming(intersect(ind,FBind)),100,'k.');
ylim([0 4]);xlim([-1,2])
hold on;line(ax,[0 1],[3 3],'Color','k');
text(ax,0.1,3.1,['pval = ',num2str(ranksum(dp_nogrooming(intersect(ind,FFind)),dp_nogrooming(intersect(ind,FBind))))])
set(ax,'XTick',0:1)
set(ax,'XTickLabel',{'FF','FB'})
%%%%%%% baseline lick latency: both FF and FB animals
% numanimals * length of timebins used for cumsum - 1 cumsum per animal
a=(cell2mat(cellfun(@(x) x(1,:),lickcumsum(ind),'UniformOutput',0)'));
figure;plot(respwindowms(1):cumsumstep:respwindowms(2)-cumsumstep, mean(a,1),'k')
hold on;patch([respwindowms(1):cumsumstep:respwindowms(2)-cumsumstep,fliplr(respwindowms(1):cumsumstep:respwindowms(2)-cumsumstep)],...
    [mean(a,1)-ci_coef*std(a)/sqrt(length(ind)),fliplr(mean(a,1)+ci_coef*std(a)/sqrt(length(ind)))],'k','EdgeAlpha',0,'FaceAlpha',0.1);
xlim([0 700])
ylabel('Lick Probability');
xlabel('time post stimulus onset (ms)');
%% laser plots: only FB animals
%%%% lick latency
FBind = find(strcmp(exptype,'FB'));
figure;
win=winter(16);
cjet=[0 0 0 ;win(end-8:end,:)];
for i = [9:-1:1]
    hold on;
    if i==3
        plot(respwindowms(1):cumsumstep:respwindowms(2)-cumsumstep,...
            mean(cell2mat(cellfun(@(x) x(i,:),lickcumsum(intersect(ind,FBind)),'UniformOutput',0)'),1),'Color',cjet(i,:),...
            'LineWidth',3)
        
    elseif i==1
        plot(respwindowms(1):cumsumstep:respwindowms(2)-cumsumstep,...
            mean(cell2mat(cellfun(@(x) x(i,:),lickcumsum(intersect(ind,FBind)),'UniformOutput',0)'),1),'Color',cjet(i,:),...
            'LineWidth',3)
    else
        plot(respwindowms(1):cumsumstep:respwindowms(2)-cumsumstep,...
            mean(cell2mat(cellfun(@(x) x(i,:),lickcumsum(intersect(ind,FBind)),'UniformOutput',0)'),1),'Color',cjet(i,:),...
            'LineWidth',1)
    end
end
ax=gca;
ax.XLabel.String = 'time (ms)';
ax.YLabel.String = 'probability of lick';

avcenters = arrayfun(@(x) num2str(round(x)),allxpoints{1},'UniformOutput',0);
avcenterslabels = cellfun(@(x) ['silencing at ',x,' ms'],avcenters,'UniformOutput',0);
legend([ax.Children(:)]',{'no silencing',avcenterslabels{:}},'FontSize',8,...
    'Location','northwest','FontWeight','normal')
legend('boxoff')
set(gcf,'Color','w');

%%% lick latency for single animals
figure;
set(gcf,'Units','normalized')
set(gcf,'Position',[0.4 0.2 0.4 0.9])

allinds = (intersect(ind,FBind)); 
for animalid_i = 1:length(allinds)
    axx=subplot(length(intersect(ind,FBind)),1,animalid_i);
    set(axx,'Units','normalized')
   
    for i = 9:-1:1
        hold on;
        if i==3
            plot(respwindowms(1):cumsumstep:respwindowms(2)-cumsumstep,...
                lickcumsum{allinds(animalid_i)}(i,:),'Color',cjet(i,:),...
                'LineWidth',2)
        elseif i==1
            plot(respwindowms(1):cumsumstep:respwindowms(2)-cumsumstep,...
                lickcumsum{allinds(animalid_i)}(i,:),'Color',cjet(i,:),...
                'LineWidth',2)
        else
            
            plot(respwindowms(1):cumsumstep:respwindowms(2)-cumsumstep,...
                lickcumsum{allinds(animalid_i)}(i,:),'Color',cjet(i,:),...
                'LineWidth',0.5)
        end
    end
    if animalid_i~=length(allinds)
        axx.XTickLabel=[];
        axx.YTickLabel=[];
    end
end
set(gcf,'Color','w');
%% %% laser plots: only FB animals: fa and miss rate
normalizeperanimal = 0;
if normalizeperanimal
    for animalid=(intersect(ind,FBind))
        allhitrate{animalid} = allhitrate{animalid} /allhitrate{animalid}(1);
    end
end

figure;
set(gcf,'Units','normalized')
set(gcf,'Position',[0.3 0.2 0.3 0.7])
axx=subplot(2,1,1);
for animalid=(intersect(ind,FBind))
    hold on
    plot(allhitrate{animalid},'.-','Color',[0.8 0.8 0.8])
end
hold on;plot(mean(cell2mat(allhitrate(intersect(ind,FBind))'),1),'k.-')
axx.YLabel.String='hit rate';
axx.XTickLabel=['no laser' ; arrayfun(@(x)  num2str(round(x)),(unique(mean(cell2mat(allxpoints')',2))),'UniformOutput',0)];
a=cell2mat(allhitrate(intersect(ind,FBind))');
[p3,~]=signrank(a(:,3),a(:,1))
if (p3>0.05);p3=nan;end
[p2,~]=signrank(a(:,2),a(:,1))
if (p2>0.05);p2=nan;end
[p4,~]=signrank(a(:,4),a(:,1));
if (p4>0.05);p4=nan;end
sigstar({[1,3],[1,2],[1,4]},[p3,p2,p4])


axx=subplot(2,1,2);
for animalid=(intersect(ind,FBind))
    hold on
    plot(allfarate{animalid},'.-','Color',[0.8 0.8 0.8])
end
hold on;plot(mean(cell2mat(allfarate(intersect(ind,FBind))'),1),'k.-')
axx.YLabel.String='false alarm rate';
axx.XTickLabel=['no laser' ; arrayfun(@(x)  num2str(round(x)),(unique(mean(cell2mat(allxpoints')',2))),'UniformOutput',0)];
set(gcf,'Color','w');
a=cell2mat(allfarate(intersect(ind,FBind))');
[p3,~]=signrank(a(:,3),a(:,1));
if (p3>0.05);p3=nan;end
[p2,~]=signrank(a(:,2),a(:,1));
if (p2>0.05);p2=nan;end
[p4,~]=signrank(a(:,4),a(:,1));
if (p4>0.05);p4=nan;end
sigstar({[1,3],[1,2],[1,4]},[p3,p2,p4])

%export_fig(gcf,sprintf('%s',['test','.png']))
%% hitrate bars
figure;
FFind = find(strcmp(exptype,'FF'));
FBind = find(strcmp(exptype,'FB'));
clear allhits
allhits(:,:,1) = cell2mat(allhitrate(intersect(ind,FFind))');
allhits(:,:,2) = cell2mat(allhitrate(intersect(ind,FBind))');
for area = 1:2
    a=allhits(:,:,area);
    subplot(1,2,area);
    b=bar(1:9,[nanmean(a(:,1)),nanmean(a(:,2)),nanmean(a(:,3)),nanmean(a(:,4)),...
        nanmean(a(:,5)),nanmean(a(:,6)),nanmean(a(:,7)),nanmean(a(:,8)),nanmean(a(:,9))]);
    [p1,~]=signrank(a(:,2),a(:,1))
    [p2,~]=signrank(a(:,3),a(:,1))
    [p3,~]=signrank(a(:,4),a(:,1))
    [p4,~]=signrank(a(:,5),a(:,1))
    b.FaceColor = 'flat';
    b.CData(1,:) = [0 0 0];
    b.CData(2,:) = [0 0 1];
    b.CData(3,:) = [0 0 1];
    b.CData(4,:) = [0 0 1];
    b.CData(5,:) = [0 0 1];
    b.CData(6,:) = [0 0 1];
    b.CData(7,:) = [0 0 1];
    b.CData(8,:) = [0 0 1];
    b.CData(9,:) = [0 0 1];
    hold on; er = errorbar(1:9,[nanmean(a(:,1)),nanmean(a(:,2)),nanmean(a(:,3)),nanmean(a(:,4)),...
        nanmean(a(:,5)),nanmean(a(:,6)),nanmean(a(:,7)),nanmean(a(:,8)),nanmean(a(:,9))],...
        [nanstd(a(:,1))/sqrt(size(a,1)), nanstd(a(:,2))/sqrt(size(a,1)),nanstd(a(:,3))/sqrt(size(a,1)),...
        nanstd(a(:,4))/sqrt(size(a,1)),nanstd(a(:,5))/sqrt(size(a,1)),nanstd(a(:,6))/sqrt(size(a,1)),...
        nanstd(a(:,7))/sqrt(size(a,1)),nanstd(a(:,8))/sqrt(size(a,1)),nanstd(a(:,9))/sqrt(size(a,1))]);
    er.Color = [1 0 0];
    er.LineStyle = 'none';
    axx=gca;
    axx.YLabel.String='hit rate';
    allXlabels = ['no laser' ; arrayfun(@(x)  num2str(round(x)),(unique(mean(cell2mat(allxpoints')',2))),'UniformOutput',0)];
    axx.XTickLabel = allXlabels(1:9);
    axx.YLim = [0 1.5];
    hold(axx,'on')
    text(0.1,1.2,sprintf('nolaser-comp-paired-pvals:\nlag1: %f\nlag2: %f\nlag3: %f\nlag4: %f',p1,p2,p3,p4))
    if area == 1
        axx.Title.String = 'feedforward';
    else
        axx.Title.String = 'feedback';
    end
end
%% distribution of lick latency for different lags - This is better substitute of slope (below section)
normalizetobs = 1;
nlags = 9;

FFind = find(strcmp(exptype,'FF'));
FBind = find(strcmp(exptype,'FB'));

alllickdelay = cell(1,nlags);
for l=1:nlags
    alllickdelay{l} = [];
end

allinds = (intersect(ind,FBind));
for animali=1:length(allinds)
    meanbs = nanmean(licklatencydist{allinds(animali)}{1});
    for l = 1:9
        if normalizetobs
            alllickdelay{l} = [alllickdelay{l} licklatencydist{allinds(animali)}{l}/meanbs];
        else
            alllickdelay{l} = [alllickdelay{l} licklatencydist{allinds(animali)}{l}];
        end
    end
end

figure;shadedErrorBar([],cellfun(@(x) nanmean(x),alllickdelay),...
    cellfun(@(x) nanstd(x)/sqrt(length(x)),alllickdelay),'k',1)
axx=gca;
axx.YLabel.String = ('mean normalized (to baseline) response latency and ste of the mean');
axx.YLabel.FontSize = 8;
axx.XTickLabel=['no laser' ; arrayfun(@(x)  num2str(round(x)),(unique(mean(cell2mat(allxpoints')',2))),'UniformOutput',0)];

figure;
for l=1:nlags
    hold on;scatter(l+rand(1,length(alllickdelay{l}))*0.2-0.1,alllickdelay{l},100,[0 0 0],'.',...
        'MarkerEdgeAlpha',0.03)
    hold on;scatter(l,nanmedian(alllickdelay{l}),100,'k','+')
end



%% lick latency at different silencing lags, quantified as the time needed to reach a certain lickrate


FFind = find(strcmp(exptype,'FF'));
FBind = find(strcmp(exptype,'FB'));

allinds = (intersect(ind,FBind));
slope = nan(length(allinds),9);
% numanimals*numlags
halftime = nan(length(allinds),9);
for animalid_i = 1:length(allinds)
   
     for i = 1:9 % lag - first one is baseline
       
        halfmahitrate = 0.5*lickcumsum{allinds(animalid_i)}(i,end);

        
        y = lickcumsum{allinds(animalid_i)}(i,:);
        x = respwindowms(1):cumsumstep:respwindowms(2)-cumsumstep;
        [~,uy_i]=unique(y);
        halftime(animalid_i,i) = interp1(y(uy_i),x(uy_i) ,halfmahitrate,'spline');
        
        % slope of the linear part
        time_e = 30;  
        time_s = 10; 
        slope(animalid_i,i) = y(time_e)-y(time_s);
     end
end
 figure;plot(slope','Color',[0.8 0.8 0.8])
 hold on;plot(nanmean(slope),'k')
 ylabel('slope (190 to 390ms)');
axx=gca;
axx.XTickLabel=['no laser' ; arrayfun(@(x)  num2str(round(x)),(unique(mean(cell2mat(allxpoints')',2))),'UniformOutput',0)];
