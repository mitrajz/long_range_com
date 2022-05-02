%%% 
% early referes to the pre decision time window: lags 1 and 2 (pre 100ms),
% pooled together
% late refers to the lags after response window onset: lags 3 to 8 (after
% 100 ms)
%
% sections 1-1 and 1-2 are for exploratory analysis. sections 2-1 and 2-2
% produce final plots (change pseudo between 0 and 1 to get onset correction or raw)
%%

% 2 crit: behavior + at least 10 trials

animallist = {'MPV17','MPV18_2',...
    'VL53','VL52','VL51','VL66',...
    'VL61','VL63','VL55','VL59',...
    'MPV33','MPV31','MPV34_2'};
preprocessinglist = {'2020_03_02_13_10_37','2020_03_02_13_27_45',...
    '2017_12_05_18_58_17','2017_12_04_15_08_47','2017_12_03_17_30_15','2018_03_20_11_53_19',...
    '2018_03_16_16_53_13','2018_04_01_12_24_38','2018_04_01_15_03_48','2018_04_01_13_24_48',...
    '2020_03_02_19_35_12','2020_03_02_19_58_50','2020_03_02_20_32_35'};
exptype = {'FB','FB',...
    'FB','FB','FB','FB',...
    'FF','FF','FF','FF',...
    'FF','FF','FF'};


allhitrate = cell(1,length(animallist));
allmissrate = cell(1,length(animallist));
allxpoints = cell(1,length(animallist));
allcrrate = cell(1,length(animallist));
allfarate = cell(1,length(animallist));
lickcumsum = cell(1,length(animallist));
licklatencydist = cell(1,length(animallist));
licklatencydist_fa = cell(1,length(animallist));

allhitn = cell(1,length(animallist));
allmissn = cell(1,length(animallist));
allcrn = cell(1,length(animallist));
allfan = cell(1,length(animallist));

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
    %%%%%%%%%%%%%%%%%%%%% response latency and hit/fa/miss/cr rates, during both baseline and laser trials
    
    min_n_trials_per_delay = 25; % 
    cumsumstep = 50; 
    onlynogrooming = 0;
    onlystableperiod = 0;
    normalizeforcorrect = 0;
    put2nanoutsidewin = 0;
    respwindowms = [100 650]; 
    onlycorrect = 1;
    
    
    [allhitrate{animal_i},allmissrate{animal_i},allfarate{animal_i},allcrrate{animal_i},allxpoints{animal_i}] = ...
        behsummaryhelper_rates(gotrialind,nogotrialind,...
        correctgotrialind,incorrectgotrialind,correctnogotrialind,incorrectnogotrialind,...
        stabletrialind,nogroomingind,LaserDelayBinned,...
        min_n_trials_per_delay,onlynogrooming,onlystableperiod,centers);
    
    [allhitn{animal_i},allmissn{animal_i},allfan{animal_i},allcrn{animal_i},allxpoints{animal_i}] = ...
        behsummaryhelper_ns(gotrialind,nogotrialind,...
        correctgotrialind,incorrectgotrialind,correctnogotrialind,incorrectnogotrialind,...
        stabletrialind,nogroomingind,LaserDelayBinned,...
        min_n_trials_per_delay,onlynogrooming,onlystableperiod,centers);
    
    
    [lickcumsum{animal_i},licklatencydist{animal_i},licklatencydist_fa{animal_i}] = behsummaryhelper_lickhist(gotrialind,nogotrialind,...
        correctgotrialind,incorrectgotrialind,correctnogotrialind,incorrectnogotrialind,...
        stabletrialind,nogroomingind,LaserDelayBinned,...
        min_n_trials_per_delay,onlynogrooming,onlystableperiod,...
        cumsumstep,normalizeforcorrect,respwindowms,put2nanoutsidewin,centers,firstlicksample,onlycorrect);
        
end
%%
 FBind =find(strcmp(exptype,'FB'))
 FFind =find(strcmp(exptype,'FF'))
 ind = 1:13;
 disp('FF animals:')
animallist(intersect(ind,FFind))
disp('FB animals:')
animallist(intersect(ind,FBind))
%% 1-1- raw plots, for inspection: d-prime comparison: pool per animals and lag 


[fb_dp,fb_delta_dp]=calc_dp(allhitn,allmissn,allcrn,allfan,intersect(ind,FBind));

figure;

subplot(1,3,1);
errorbar(1:3,[mean(fb_dp.bs),mean(fb_dp.early),mean(fb_dp.late)],[std(fb_dp.bs)/sqrt(numel(fb_dp.bs)),std(fb_dp.early)/sqrt(numel(fb_dp.early)),std(fb_dp.late)/sqrt(numel(fb_dp.late))])
ylim([0 3])
xlim([0 4])
ylabel('d prime')
text(2,2.9,['ranksum early to bs, p = ',num2str(ranksum(fb_dp.bs,fb_dp.early))],'HorizontalAlignment','center');
text(3,2.7,['ranksum late to bs, p = ',num2str(ranksum(fb_dp.bs,fb_dp.late))],'HorizontalAlignment','center');

subplot(1,3,2);
hold on;bar(1,mean(fb_delta_dp.early));
hold on;bar(2,mean(fb_delta_dp.late));
hold on;errorbar(1:2,[mean(fb_delta_dp.early),mean(fb_delta_dp.late)],[std(fb_delta_dp.early)/sqrt(numel(fb_delta_dp.early)),std(fb_delta_dp.late)/sqrt(numel(fb_delta_dp.late))])
ylim([-1 1])
xlim([0 3])
ylabel('delta d prime')
text(1,0.9,['signrank early to 0 = ',num2str(signrank(fb_delta_dp.early)) ',n=',num2str(numel(fb_delta_dp.early)) ]);
text(2,0.8,['signrank late to 0 = ',num2str(signrank(fb_delta_dp.late)) ',n=',num2str(numel(fb_delta_dp.late)) ]);

f=gcf;f.Children(1).Title.String = 'fb';

subplot(1,3,3);
hold on;bar(0.9,mean(repmat(fb_dp.bs,2,1)),0.2,'FaceColor','k');hold on;bar(1.1,mean(fb_dp.early),0.2,'FaceColor','b')
hold on;bar(1.9,mean(repmat(fb_dp.bs,6,1)),0.2,'FaceColor','k');hold on;bar(2.1,mean(fb_dp.late),0.2,'FaceColor','b')
hold on;errorbar([0.9,1.1],[mean(repmat(fb_dp.bs,2,1)),mean(fb_dp.early)],...
    [std(repmat(fb_dp.bs,2,1))/sqrt(numel(repmat(fb_dp.bs,2,1))),std(fb_dp.early)/sqrt(numel(fb_dp.early))])
hold on;errorbar([1.9,2.1],[mean(repmat(fb_dp.bs,6,1)),mean(fb_dp.late)],...
    [std(repmat(fb_dp.bs,6,1))/sqrt(numel(repmat(fb_dp.bs,6,1))),std(fb_dp.late)/sqrt(numel(fb_dp.late))])
ylim([0 3])
xlim([0 3])
ylabel('d prime')
text(1,2.5,['signrank early to 0 = ',num2str(signrank(repmat(fb_dp.bs,2,1),fb_dp.early)), ',n=',num2str(numel(fb_dp.early)) ]);
text(2,2.7,['signrank late to 0 = ',num2str(signrank(repmat(fb_dp.bs,6,1),fb_dp.late)) ',n=',num2str(numel(fb_dp.late)) ]);

%%%%%%%%%% FF: 
[ff_dp,ff_delta_dp]=calc_dp(allhitn,allmissn,allcrn,allfan,intersect(ind,FFind));

figure;

subplot(1,3,1);errorbar(1:3,[mean(ff_dp.bs),mean(ff_dp.early),mean(ff_dp.late)],[std(ff_dp.bs)/sqrt(numel(ff_dp.bs)),std(ff_dp.early)/sqrt(numel(ff_dp.early)),std(ff_dp.late)/sqrt(numel(ff_dp.late))])
ylim([0 3])
xlim([0 4])
ylabel('d prime')
text(2,2.9,['ranksum early to bs, p = ',num2str(ranksum(ff_dp.bs,ff_dp.early))],'HorizontalAlignment','center');
text(3,2.7,['ranksum late to bs, p = ',num2str(ranksum(ff_dp.bs,ff_dp.late))],'HorizontalAlignment','center');

subplot(1,3,2);
hold on;bar(1,mean(ff_delta_dp.early));
hold on;bar(2,mean(ff_delta_dp.late));
hold on;errorbar(1:2,[mean(ff_delta_dp.early),mean(ff_delta_dp.late)],[std(ff_delta_dp.early)/sqrt(numel(ff_delta_dp.early)),std(ff_delta_dp.late)/sqrt(numel(ff_delta_dp.late))])
hold on;scatter(ones(1,length(ff_delta_dp.early)),ff_delta_dp.early)
hold on;scatter(1+ones(1,length(ff_delta_dp.late)),ff_delta_dp.late)

ylim([-1 1])
xlim([0 3])
ylabel('delta d prime')
text(1,0.9,['signrank early to 0 = ',num2str(signrank(ff_delta_dp.early)), ',n=',num2str(numel(ff_delta_dp.early)) ]);
text(2,0.8,['signrank late to 0 = ',num2str(signrank(ff_delta_dp.late)) ',n=',num2str(numel(ff_delta_dp.late)) ]);

f=gcf;f.Children(1).Title.String = 'ff';

% final one
subplot(1,3,3);
hold on;bar(0.9,mean(repmat(ff_dp.bs,2,1)),0.2,'FaceColor','k');hold on;bar(1.1,mean(ff_dp.early),0.2,'FaceColor','b')
hold on;bar(1.9,mean(repmat(ff_dp.bs,6,1)),0.2,'FaceColor','k');hold on;bar(2.1,mean(ff_dp.late),0.2,'FaceColor','b')
hold on;errorbar([0.9,1.1],[mean(repmat(ff_dp.bs,2,1)),mean(ff_dp.early)],...
    [std(repmat(ff_dp.bs,2,1))/sqrt(numel(repmat(ff_dp.bs,2,1))),std(ff_dp.early)/sqrt(numel(ff_dp.early))])
hold on;errorbar([1.9,2.1],[mean(repmat(ff_dp.bs,6,1)),mean(ff_dp.late)],...
    [std(repmat(ff_dp.bs,6,1))/sqrt(numel(repmat(ff_dp.bs,6,1))),std(ff_dp.late)/sqrt(numel(ff_dp.late))])
ylim([0 3])
xlim([0 3])
ylabel('d prime')
text(1,2,['signrank early to 0 = ',num2str(signrank(repmat(ff_dp.bs,2,1),ff_dp.early)), ',n=',num2str(numel(ff_dp.early)) ]);
text(2,2.2,['signrank late to 0 = ',num2str(signrank(repmat(ff_dp.bs,6,1),ff_dp.late)) ',n=',num2str(numel(ff_dp.late)) ]);

%% 1-2- raw plots, for inspection: reaction time

%%% toggle exptype between ff and fb to plots both

exptype = 'FB';

if strcmp(exptype,'FF')
    EXPind = intersect(ind,FFind)
elseif strcmp(exptype,'FB')
    EXPind = intersect(ind,FBind)
else
    error('specify exptype')
end

normalizetobs = 1; % this normalizes per animal to take care of per animal differences in average response time

alllickdelay = cell(1,9);
for l=1:9
    alllickdelay{l} = [];
end

for animali=1:length(EXPind)
    meanbs = nanmean(licklatencydist{EXPind(animali)}{1});
    for l = 1:9
        if normalizetobs
            alllickdelay{l} = [alllickdelay{l} 100*(licklatencydist{EXPind(animali)}{l}-meanbs)/meanbs];
        else
            alllickdelay{l} = [alllickdelay{l} licklatencydist{EXPind(animali)}{l}];
        end
    end
end
% mean and standard error of the mean - average normalized (to baseline) latency 

figure;subplot(1,3,1);shadedErrorBar([],cellfun(@(x) nanmean(x),alllickdelay),...
    cellfun(@(x) nanstd(x)/sqrt(length(x)),alllickdelay),'k',1)
axx=gca;
axx.YLabel.String = ('mean normalized (to baseline) response latency and ste of the mean');
axx.YLabel.FontSize = 8;
axx.XTickLabel=['no laser' ; arrayfun(@(x)  num2str(round(x)),(unique(mean(cell2mat(allxpoints')',2))),'UniformOutput',0)];



clear rt
rt.bs = alllickdelay{1};
rt.early = [alllickdelay{2},alllickdelay{3}];
rt.late = [alllickdelay{4},alllickdelay{5},alllickdelay{6},alllickdelay{7},alllickdelay{8},alllickdelay{9}];
subplot(1,3,2);
errorbar(1:3,[mean(rt.bs),mean(rt.early),mean(rt.late)],[std(rt.bs)/sqrt(numel(rt.bs)),std(rt.early)/sqrt(numel(rt.early)),std(rt.late)/sqrt(numel(rt.late))])
xlim([0 4])

text(2,0,['ranksum early to bs = ',num2str(ranksum(rt.bs,rt.early))],'HorizontalAlignment','center')
text(2,1,['signrank early to 0 = ',num2str(signrank(rt.early)),',n=',num2str(numel(rt.early))],'HorizontalAlignment','center')
text(3,2,['ranksum late to bs = ',num2str(ranksum(rt.bs,rt.late))],'HorizontalAlignment','center')
text(3,3,['signrank late to 0 = ',num2str(signrank(rt.late)),',n=',num2str(numel(rt.late))],'HorizontalAlignment','center')

subplot(1,3,3);
histogram(rt.bs,-500:42:500,'FaceColor','k','FaceAlpha',0.2,'Normalization','pdf')
hold on;histogram(rt.early,-500:42:500,'FaceColor','b','FaceAlpha',0.2,'Normalization','pdf')
hold on;histogram(rt.late,-500:42:500,'FaceColor','r','FaceAlpha',0.2,'Normalization','pdf')
%% 2-1 - processed plots, for reporting: d prime but with lick onset after bin start cntrl
%%% This section produces %change in d prime (new-bs)/bs where bs is
%%% baseline d prime per animal. psuedo gives the option to either only
%%% count licks after the bin onset, or all of them. the 2 columns are
%%% baseline and laser, reported for early and late separately

%%% toggle exptype between ff and fb to plots both

exptype = 'FB';
psuedo = 1;


if strcmp(exptype,'FF')
    EXPind = intersect(ind,FFind)
elseif strcmp(exptype,'FB')
    EXPind = intersect(ind,FBind)
else
    error('specify exptype')
end

bs_go_lickafterlaser_rate = nan(numel(EXPind),8);
ls_go_lickafterlaser_rate = nan(numel(EXPind),8);
bs_nogo_lickafterlaser_rate = nan(numel(EXPind),8);
ls_nogo_lickafterlaser_rate = nan(numel(EXPind),8);
    
for i = 1:numel(EXPind)
    animali = EXPind(i);
    % go

    bs_go_n = allhitn{animali}(1) + allmissn{animali}(1);
    for l = 1:8       
        % bs corresponding to laser
        if psuedo
            bs_go_lickafterlaser = numel(find(licklatencydist{animali}{1}>allxpoints{animali}(l)));
        else
            bs_go_lickafterlaser = numel(licklatencydist{animali}{1});
        end
        if bs_go_lickafterlaser/bs_go_n == 0
            bs_go_lickafterlaser_rate(i,l) = 0.5/(bs_go_n-bs_go_lickafterlaser); % the denum is the equivalent to the number of miss trials
        elseif bs_go_lickafterlaser/bs_go_n == 1
            bs_go_lickafterlaser_rate(i,l) = (bs_go_lickafterlaser-0.5)/bs_go_lickafterlaser;
        else
            bs_go_lickafterlaser_rate(i,l) = bs_go_lickafterlaser/bs_go_n; 
        end
        
        % laser       
        ls_go_n= allhitn{animali}(l+1) + allmissn{animali}(l+1);
        if psuedo
            ls_go_lickafterlaser = numel(find(licklatencydist{animali}{l+1}>allxpoints{animali}(l)));
        else
            ls_go_lickafterlaser = numel(licklatencydist{animali}{l+1});
        end
         if ls_go_lickafterlaser/ls_go_n == 0
            ls_go_lickafterlaser_rate(i,l) = 0.5/(ls_go_n-ls_go_lickafterlaser); % the denum is the equivalent to the number of miss trials
        elseif ls_go_lickafterlaser/ls_go_n == 1
            ls_go_lickafterlaser_rate(i,l) = (ls_go_lickafterlaser-0.5)/ls_go_lickafterlaser;
        else
            ls_go_lickafterlaser_rate(i,l) = ls_go_lickafterlaser/ls_go_n; 
         end
    end
    
    % nogo

    bs_nogo_n = allcrn{animali}(1) + allfan{animali}(1);
    for l = 1:8
       
        % bs corresponding to laser
        if psuedo
            bs_nogo_lickafterlaser = numel(find(licklatencydist_fa{animali}{1}>allxpoints{animali}(l)));
        else
            bs_nogo_lickafterlaser = numel(licklatencydist_fa{animali}{1});
        end
        if bs_nogo_lickafterlaser/bs_nogo_n == 0
            bs_nogo_lickafterlaser_rate(i,l) = 0.5/(bs_nogo_n-bs_nogo_lickafterlaser); % the denum is the equivalent to the number of cr trials
        elseif bs_nogo_lickafterlaser/bs_nogo_n == 1
            bs_nogo_lickafterlaser_rate(i,l) = (bs_nogo_lickafterlaser-0.5)/bs_nogo_lickafterlaser;
        else
            bs_nogo_lickafterlaser_rate(i,l) = bs_nogo_lickafterlaser/bs_nogo_n; 
        end
    
        
        % laser
        ls_nogo_n= allcrn{animali}(l+1) + allfan{animali}(l+1);
        if psuedo
            ls_nogo_lickafterlaser = numel(find(licklatencydist_fa{animali}{l+1}>allxpoints{animali}(l)));
        else
            ls_nogo_lickafterlaser = numel(licklatencydist_fa{animali}{l+1});
        end
         if ls_nogo_lickafterlaser/ls_nogo_n == 0
            ls_nogo_lickafterlaser_rate(i,l) = 0.5/(ls_nogo_n-ls_nogo_lickafterlaser); % the denum is the equivalent to the number of cr trials
        elseif ls_nogo_lickafterlaser/ls_nogo_n == 1
            ls_nogo_lickafterlaser_rate(i,l) = (ls_nogo_lickafterlaser-0.5)/ls_nogo_lickafterlaser;
         else
            ls_nogo_lickafterlaser_rate(i,l) = ls_nogo_lickafterlaser/ls_nogo_n; 
        end
    end
    
end
bs_dprime_perbin = norminv(bs_go_lickafterlaser_rate) - norminv(bs_nogo_lickafterlaser_rate);
ls_dprime_perbin = norminv(ls_go_lickafterlaser_rate) - norminv(ls_nogo_lickafterlaser_rate);


ls_dprime_perbin = 100*(ls_dprime_perbin - repmat(bs_dprime_perbin(:,1),1,8))./repmat(bs_dprime_perbin(:,1),1,8);
bs_dprime_perbin = 100*(bs_dprime_perbin - repmat(bs_dprime_perbin(:,1),1,8))./repmat(bs_dprime_perbin(:,1),1,8);



bs_dprime_perbin_early = reshape(bs_dprime_perbin(:,1:2),1,[]);
bs_dprime_perbin_late = reshape(bs_dprime_perbin(:,3:8),1,[]);
ls_dprime_perbin_early = reshape(ls_dprime_perbin(:,1:2),1,[]);
ls_dprime_perbin_late = reshape(ls_dprime_perbin(:,3:8),1,[]);

figure;
hold on;bar(0.9,mean(bs_dprime_perbin_early),0.2,'FaceColor','k','FaceAlpha',0.5);
hold on;bar(1.1,mean(ls_dprime_perbin_early),0.2,'FaceColor','b','FaceAlpha',0.5);
hold on;line([0.9,1.1],[mean(bs_dprime_perbin_early);mean(ls_dprime_perbin_early)]);

hold on;bar(1.9,mean(bs_dprime_perbin_late),0.2,'FaceColor','k','FaceAlpha',0.5);
hold on;bar(2.1,mean(ls_dprime_perbin_late),0.2,'FaceColor','b','FaceAlpha',0.5);
hold on;line([1.9,2.1],[mean(bs_dprime_perbin_late);mean(ls_dprime_perbin_late)]);

hold on;errorbar([0.9,1.9],[mean(bs_dprime_perbin_early),mean(bs_dprime_perbin_late)],...
    [std(bs_dprime_perbin_early)/sqrt(numel(bs_dprime_perbin_early)),std(bs_dprime_perbin_late)/sqrt(numel(bs_dprime_perbin_late))],'k')
hold on;errorbar([1.1,2.1],[mean(ls_dprime_perbin_early),mean(ls_dprime_perbin_late)],...
    [std(ls_dprime_perbin_early)/sqrt(numel(ls_dprime_perbin_early)),std(ls_dprime_perbin_late)/sqrt(numel(ls_dprime_perbin_late))],'b')
xlim([0 3])
ylim([-70 3])
text(1,2.9,['signrank bs ls early = ',num2str(signrank(bs_dprime_perbin_early,ls_dprime_perbin_early)),',n=',num2str(numel(bs_dprime_perbin_early))]);
text(2,2.8,['signrank bs ls late = ',num2str(signrank(bs_dprime_perbin_late,ls_dprime_perbin_late)),',n=',num2str(numel(bs_dprime_perbin_late))]);
title('')
ylabel(['performance change (d-prime or pseudo d-prime) - ',exptype])
 

%% 2-2 - processed plots, for reporting: reaction time, but with lick onset after bin start cntrl

%%% toggle exptype between ff and fb to plots both

exptype = 'FB';
psuedo = 0;
extraplots = 0;

if strcmp(exptype,'FF')
    EXPind = intersect(ind,FFind)
elseif strcmp(exptype,'FB')
    EXPind = intersect(ind,FBind)
else
    error('specify exptype')
end

alllickdelay = cell(1,9);
for l=1:9
    alllickdelay{l} = [];
end


alllickdelay_bs = cell(1,9);
for l=1:9
    alllickdelay_bs{l} = [];
end

for animali=1:length(EXPind)
    meanbs = nanmean(licklatencydist{EXPind(animali)}{1});
    for l = 1:8
        if psuedo
            afterlaser = licklatencydist{EXPind(animali)}{l+1}(find(licklatencydist{EXPind(animali)}{l+1}>allxpoints{EXPind(animali)}(l)));
            afterlaser_bs = licklatencydist{EXPind(animali)}{1}(find(licklatencydist{EXPind(animali)}{1}>allxpoints{EXPind(animali)}(l)));
        else
            afterlaser = licklatencydist{EXPind(animali)}{l+1};
            afterlaser_bs = licklatencydist{EXPind(animali)}{1};
        end
         alllickdelay{l+1} = [alllickdelay{l+1} 100*(afterlaser-meanbs)/meanbs];
         alllickdelay_bs{l+1} = [alllickdelay_bs{l+1} 100*(afterlaser_bs-meanbs)/meanbs];
     
    end
end
alllickdelay(1)=[];
alllickdelay_bs(1)=[];

clear ls_rt bs_rt

ls_rt.early = [alllickdelay{1},alllickdelay{2}];
ls_rt.late = [alllickdelay{3},alllickdelay{4},alllickdelay{5},alllickdelay{6},alllickdelay{7},alllickdelay{8}];

bs_rt.early = [alllickdelay_bs{1},alllickdelay_bs{2}];
bs_rt.late = [alllickdelay_bs{3},alllickdelay_bs{4},alllickdelay_bs{5},alllickdelay_bs{6},alllickdelay_bs{7},alllickdelay_bs{8}];

figure;
subplot(1,3,2);
hold on; bar(0.9,mean(bs_rt.early),0.2,'FaceColor','k','FaceAlpha',0.5);
hold on; bar(1.1,mean(ls_rt.early),0.2,'FaceColor','b','FaceAlpha',0.5);
hold on; bar(1.9,mean(bs_rt.late),0.2,'FaceColor','k','FaceAlpha',0.5);
hold on; bar(2.1,mean(ls_rt.late),0.2,'FaceColor','b','FaceAlpha',0.5);
hold on;errorbar([0.9,1.9],[mean(bs_rt.early),mean(bs_rt.late)],[std(bs_rt.early)/sqrt(numel(bs_rt.early)),std(bs_rt.late)/sqrt(numel(bs_rt.late))],'k')
hold on;errorbar([1.1,2.1],[mean(ls_rt.early),mean(ls_rt.late)],[std(ls_rt.early)/sqrt(numel(ls_rt.early)),std(ls_rt.late)/sqrt(numel(ls_rt.late))],'b')
xlim([0 4])
%%%% 
text(1,0,['ranksum early ls to bs = ',num2str(ranksum(bs_rt.early,ls_rt.early))],'HorizontalAlignment','center')
text(2,2,['ranksum late ls to bs = ',num2str(ranksum(bs_rt.late ,ls_rt.late ))],'HorizontalAlignment','center')
ylim([-5 20])
title(exptype)


if extraplots
    subplot(1,3,1);shadedErrorBar([],cellfun(@(x) nanmean(x),alllickdelay_bs),...
        cellfun(@(x) nanstd(x)/sqrt(length(x)),alllickdelay_bs),'k',1)
    hold on;subplot(1,3,1);shadedErrorBar([],cellfun(@(x) nanmean(x),alllickdelay),...
        cellfun(@(x) nanstd(x)/sqrt(length(x)),alllickdelay),'b',1)
    axx=gca;
  
    axx.YLabel.FontSize = 8;
    axx.XTick = 1:8;
    axx.XLim = [0,9];
    axx.XTickLabel=[arrayfun(@(x)  num2str(round(x)),(unique(mean(cell2mat(allxpoints')',2))),'UniformOutput',0)];
    
    subplot(1,3,3); 
    hold on;histogram(bs_rt.early,-500:42:500,'FaceColor','k','FaceAlpha',0.2,'Normalization','pdf')
    hold on;histogram(ls_rt.early,-500:42:500,'FaceColor','b','FaceAlpha',0.2,'Normalization','pdf')
    
end

%% helper functions
function [dp,delta_dp]=calc_dp(allhitn,allmissn,allcrn,allfan,indx)
% dim:nanimals*nlags+1
allhitmat = cell2mat(allhitn(indx)');
allmissmat = cell2mat(allmissn(indx)');
allfamat = cell2mat(allfan(indx)');
allcrmat = cell2mat(allcrn(indx)');


% delta dp is difference wrt the average bs condition
hr_perbin = allhitmat./(allhitmat+allmissmat);
fa_perbin = allfamat./(allfamat+allcrmat);

% fixes 0s and 1s for hit ratre and fa rate
hr_perbin(find(hr_perbin==1)) = (allhitmat(find(hr_perbin==1))  - 0.5)./allhitmat(find(hr_perbin==1));
hr_perbin(find(hr_perbin==0)) = (0.5)./allmissmat(find(hr_perbin==0)); 

fa_perbin(find(fa_perbin==1)) = (allfamat(find(fa_perbin==1))  - 0.5)./allfamat(find(fa_perbin==1));
fa_perbin(find(fa_perbin==0)) = (0.5)./allcrmat(find(fa_perbin==0)); 


dprime_perbin = norminv(hr_perbin) - norminv(fa_perbin);

dp.bs = dprime_perbin(:,1);
dp.early = reshape(dprime_perbin(:,2:3),1,[]);
dp.late = reshape(dprime_perbin(:,4:9),1,[]);


delta = dprime_perbin - repmat(dprime_perbin(:,1),1,9) ;

delta_dp.early = reshape(delta(:,2:3),1,[]);
delta_dp.late = reshape(delta(:,4:9),1,[]);


end

