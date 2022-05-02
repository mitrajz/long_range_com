%% first get the resting firing rate of cells
bsfr = nan(1,length(V1cells));
for i=1:length(V1cells)
    bsfr(i) = ...
        nanmean([V1cells{i}.nbsAv.go(1:floor(length(V1cells{i}.nbsAv.go)/2)),V1cells{i}.nbsAv.nogo(1:floor(length(V1cells{i}.nbsAv.nogo)/2))]);
end
%% 1. it is not causal yet (solved), 2.arbitrary threshold of firing 3.sem?
figure;
for lag=1:8
    V1sel = nan(length(V1cells),length(V1cells{1}.nbsAv.go));
    V1selslc = nan(length(V1cells),length(V1cells{1}.nbsAv.go));
    clear feat feat2
    causalbinning=1;
    
    for i=1:length(V1cells)
        highratetimeind = (V1cells{i}.nbsAv.go+V1cells{i}.nbsAv.nogo)>0.25;
        %V1sel(i,:) = ((V1cells{i}.nbsAv.go-bsfr(i))/bsfr(i) - (V1cells{i}.nbsAv.nogo-bsfr(i))/bsfr(i))./((V1cells{i}.nbsAv.go-bsfr(i))/bsfr(i) + (V1cells{i}.nbsAv.nogo-bsfr(i))/bsfr(i));      V1sel(i,:) = ((V1cells{i}.nbsAv.go-bsfr(i))/bsfr(i) - (V1cells{i}.nbsAv.nogo-bsfr(i))/bsfr(i))./((V1cells{i}.nbsAv.go-bsfr(i))/bsfr(i) + (V1cells{i}.nbsAv.nogo-bsfr(i))/bsfr(i));
        V1sel(i,:) = (V1cells{i}.nbsAv.go - V1cells{i}.nbsAv.nogo)./(V1cells{i}.nbsAv.go + V1cells{i}.nbsAv.nogo);
        
        V1sel(i,find(highratetimeind == 0)) = nan;
        
        highratetimeind = (V1cells{i}.nlsAv.go(lag,:)+V1cells{i}.nlsAv.nogo(lag,:))>0.25;
        V1selslc(i,:) = (V1cells{i}.nlsAv.go(lag,:)- V1cells{i}.nlsAv.nogo(lag,:))./((V1cells{i}.nlsAv.go(lag,:) + V1cells{i}.nlsAv.nogo(lag,:)));
        V1selslc(i,find(highratetimeind == 0)) = nan;
        
    end
    
    if causalbinning
        V1sel = circshift(V1sel,1,2);
        V1selslc = circshift(V1selslc,1,2);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(8,1,lag);t1=shadedErrorBar(((-Window):1*edgestep:(Window-1))/(samplingf/1000),nanmean(V1sel),1.96*nanstd(V1sel)/sqrt(length(V1cells)),{'LineWidth',1.6,'Color',[0 0 0]},1);
    hold on;t2=shadedErrorBar(((-Window):1*edgestep:(Window-1))/(samplingf/1000),nanmean(V1selslc),1.96*nanstd(V1selslc)/sqrt(length(V1cells)),{'LineWidth',1.6,'Color',[0 0 1]},1);
    xlim([-400 1000]);ylim([-0.3 0.3])
    axx=gca;
    hold on; patch([0 500 500 0],[axx.YLim(1) axx.YLim(1) axx.YLim(2) axx.YLim(2)],[1 0.7 0.7],...
        'FaceAlpha',0.1,'EdgeColor','none')
    
    startpoint=V1cells{1}.smb_centers.go(lag);
    dur=40;
    dur2=75;
    dur3=150;
    % hold on;patch([ startpoint startpoint+dur startpoint+dur startpoint],[0 0 10 10],'b','FaceAlpha',0.1)
    hold on;patch([ startpoint startpoint+dur3 startpoint+dur3 startpoint],[axx.YLim(1) axx.YLim(1) axx.YLim(2) axx.YLim(2)],[0 0.7 1],'FaceAlpha',0.1,'EdgeAlpha',0.5)
    hold on;patch([ startpoint startpoint+dur2 startpoint+dur2 startpoint],[axx.YLim(1) axx.YLim(1) axx.YLim(2) axx.YLim(2)],[0 0.7 1],'FaceAlpha',0.1,'EdgeAlpha',0)
    
    % lgd = legend([t1.mainLine t2.mainLine],'selectivity index, baseline','selectivity index, with silencing',...
    %     'Location','southeast');
    % lgd.FontSize = 7;
    if lag ==8
        xlabel('time (ms)');
        lgd = legend([t1.mainLine t2.mainLine],'selectivity index, baseline','selectivity index, with silencing',...
            'Location','northwest');
        lgd.FontSize = 7;
    else
        axx.XTickLabel = [];
        axx.YTickLabel = [];
    end
    set(axx,'Position',[0.13,(0.95-0.1*lag),0.775,0.1])
end

%% selectivity index for all cells vs.indgo
lag=2
V1sel = nan(length(V1cells),length(V1cells{1}.nbsAv.go));
V1selslc = nan(length(V1cells),length(V1cells{1}.nbsAv.go));
clear feat feat2
causalbinning=1;

for i=1:length(V1cells)
    highratetimeind = (V1cells{i}.nbsAv.go+V1cells{i}.nbsAv.nogo)>0.25;
    %V1sel(i,:) = ((V1cells{i}.nbsAv.go-bsfr(i))/bsfr(i) - (V1cells{i}.nbsAv.nogo-bsfr(i))/bsfr(i))./((V1cells{i}.nbsAv.go-bsfr(i))/bsfr(i) + (V1cells{i}.nbsAv.nogo-bsfr(i))/bsfr(i));      V1sel(i,:) = ((V1cells{i}.nbsAv.go-bsfr(i))/bsfr(i) - (V1cells{i}.nbsAv.nogo-bsfr(i))/bsfr(i))./((V1cells{i}.nbsAv.go-bsfr(i))/bsfr(i) + (V1cells{i}.nbsAv.nogo-bsfr(i))/bsfr(i));
    V1sel(i,:) = (V1cells{i}.nbsAv.go - V1cells{i}.nbsAv.nogo)./(V1cells{i}.nbsAv.go + V1cells{i}.nbsAv.nogo);
    
    V1sel(i,find(highratetimeind == 0)) = nan;
    
    highratetimeind = (V1cells{i}.nlsAv.go(lag,:)+V1cells{i}.nlsAv.nogo(lag,:))>0.25;
    V1selslc(i,:) = (V1cells{i}.nlsAv.go(lag,:)- V1cells{i}.nlsAv.nogo(lag,:))./((V1cells{i}.nlsAv.go(lag,:) + V1cells{i}.nlsAv.nogo(lag,:)));
    V1selslc(i,find(highratetimeind == 0)) = nan;
    
end

if causalbinning
    V1sel = circshift(V1sel,1,2);
    V1selslc = circshift(V1selslc,1,2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
subplot(2,1,1);t1=shadedErrorBar(((-Window):1*edgestep:(Window-1))/(samplingf/1000),nanmean(V1sel),1.96*nanstd(V1sel)/sqrt(length(V1cells)),{'LineWidth',1.6,'Color',[0 0 0]},1);
hold on;t2=shadedErrorBar(((-Window):1*edgestep:(Window-1))/(samplingf/1000),nanmean(V1selslc),1.96*nanstd(V1selslc)/sqrt(length(V1cells)),{'LineWidth',1.6,'Color',[0 0 1]},1);
xlim([-400 1000]);ylim([-0.5 0.5])
axx=gca;
hold on; patch([0 500 500 0],[axx.YLim(1) axx.YLim(1) axx.YLim(2) axx.YLim(2)],[1 0.7 0.7],...
    'FaceAlpha',0.1,'EdgeColor','none')

startpoint=V1cells{1}.smb_centers.go(lag);
dur=40;
dur2=75;
dur3=150;
% hold on;patch([ startpoint startpoint+dur startpoint+dur startpoint],[0 0 10 10],'b','FaceAlpha',0.1)
hold on;patch([ startpoint startpoint+dur3 startpoint+dur3 startpoint],[axx.YLim(1) axx.YLim(1) axx.YLim(2) axx.YLim(2)],[0 0.7 1],'FaceAlpha',0.1,'EdgeAlpha',0.5)
hold on;patch([ startpoint startpoint+dur2 startpoint+dur2 startpoint],[axx.YLim(1) axx.YLim(1) axx.YLim(2) axx.YLim(2)],[0 0.7 1],'FaceAlpha',0.1,'EdgeAlpha',0)

lgd = legend([t1.mainLine t2.mainLine],'selectivity index, baseline','selectivity index, with silencing',...
    'Location','southeast');
lgd.FontSize = 7;
xlabel('time (ms)');

set(axx,'Position',[0.13,0.5,0.8,0.4]);
axx.XTickLabel = [];
axx.YTickLabel = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V1selslc = V1selslc(indgo_g_40,:);
V1sel = V1sel(indgo_g_40,:);

subplot(2,1,2);t1=shadedErrorBar(((-Window):1*edgestep:(Window-1))/(samplingf/1000),nanmean(V1sel),1.96*nanstd(V1sel)/sqrt(length(V1cells)),{'LineWidth',1.6,'Color',[0 0 0]},1);
hold on;t2=shadedErrorBar(((-Window):1*edgestep:(Window-1))/(samplingf/1000),nanmean(V1selslc),1.96*nanstd(V1selslc)/sqrt(length(V1cells)),{'LineWidth',1.6,'Color',[0 0 1]},1);
xlim([-400 1000]);ylim([-0.5 0.5])
axx=gca;
hold on; patch([0 500 500 0],[axx.YLim(1) axx.YLim(1) axx.YLim(2) axx.YLim(2)],[1 0.7 0.7],...
    'FaceAlpha',0.1,'EdgeColor','none')

startpoint=V1cells{1}.smb_centers.go(lag);
dur=40;
dur2=75;
dur3=150;
% hold on;patch([ startpoint startpoint+dur startpoint+dur startpoint],[0 0 10 10],'b','FaceAlpha',0.1)
hold on;patch([ startpoint startpoint+dur3 startpoint+dur3 startpoint],[axx.YLim(1) axx.YLim(1) axx.YLim(2) axx.YLim(2)],[0 0.7 1],'FaceAlpha',0.1,'EdgeAlpha',0.5)
hold on;patch([ startpoint startpoint+dur2 startpoint+dur2 startpoint],[axx.YLim(1) axx.YLim(1) axx.YLim(2) axx.YLim(2)],[0 0.7 1],'FaceAlpha',0.1,'EdgeAlpha',0)

lgd = legend([t1.mainLine t2.mainLine],'selectivity index, baseline','selectivity index, with silencing',...
    'Location','southeast');
lgd.FontSize = 7;
xlabel('time (ms)');

set(axx,'Position',[0.13,0.1,0.8,0.4])
%% selectivity baseline vs laser : this one is not causal, evaluated at the through point

osi1 = nan(size(V1cells));
osi2 = nan(size(V1cells));

lag=2;
tp=28;
for i = 1:length(V1cells)
    highratetimeind = (V1cells{i}.nbsAv.go+V1cells{i}.nbsAv.nogo)>0.25;
    if highratetimeind(tp)>0
        osi1(i) = (V1cells{i}.nbsAv.go(tp) - V1cells{i}.nbsAv.nogo(tp))./(V1cells{i}.nbsAv.go(tp) + V1cells{i}.nbsAv.nogo(tp));
    end
    
    highratetimeind = (V1cells{i}.nlsAv.go(lag,:)+V1cells{i}.nlsAv.nogo(lag,:))>0.25;
    if highratetimeind(tp)>0
        osi2(i) = (V1cells{i}.nlsAv.go(lag,tp)- V1cells{i}.nlsAv.nogo(lag,tp))./((V1cells{i}.nlsAv.go(lag,tp) + V1cells{i}.nlsAv.nogo(lag,tp)));
    end
end

figure;scatter(osi1,osi2,200,'.');hold on;line([-2 2],[-2 2])
xlabel('selectivity index baseline');ylabel('selectivity index with silencing');
%% grooming vs. nogrooming 

lag=2
V1sel = nan(length(V1cells),length(V1cells{1}.nbsAv.go));
V1selslc = nan(length(V1cells),length(V1cells{1}.nbsAv.go));
clear feat feat2
causalbinning=1;

for i=1:length(V1cells)
    highratetimeind = (V1cells{i}.nbsAv.go+V1cells{i}.nbsAv.nogo)>0.25;
    %V1sel(i,:) = ((V1cells{i}.nbsAv.go-bsfr(i))/bsfr(i) - (V1cells{i}.nbsAv.nogo-bsfr(i))/bsfr(i))./((V1cells{i}.nbsAv.go-bsfr(i))/bsfr(i) + (V1cells{i}.nbsAv.nogo-bsfr(i))/bsfr(i));      V1sel(i,:) = ((V1cells{i}.nbsAv.go-bsfr(i))/bsfr(i) - (V1cells{i}.nbsAv.nogo-bsfr(i))/bsfr(i))./((V1cells{i}.nbsAv.go-bsfr(i))/bsfr(i) + (V1cells{i}.nbsAv.nogo-bsfr(i))/bsfr(i));
    V1sel(i,:) = (V1cells{i}.nbsAv.go - V1cells{i}.nbsAv.nogo)./(V1cells{i}.nbsAv.go + V1cells{i}.nbsAv.nogo);
    
    V1sel(i,find(highratetimeind == 0)) = nan;
    
    highratetimeind = (V1cells{i}.nlsAv.go(lag,:)+V1cells{i}.nlsAv.nogo(lag,:))>0.25;
    V1selslc(i,:) = (V1cells{i}.nlsAv.go(lag,:)- V1cells{i}.nlsAv.nogo(lag,:))./((V1cells{i}.nlsAv.go(lag,:) + V1cells{i}.nlsAv.nogo(lag,:)));
    V1selslc(i,find(highratetimeind == 0)) = nan;
    
end

if causalbinning
    V1sel = circshift(V1sel,1,2);
    V1selslc = circshift(V1selslc,1,2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
t1=shadedErrorBar(((-Window):1*edgestep:(Window-1))/(samplingf/1000),nanmean(V1sel),1.96*nanstd(V1sel)/sqrt(length(V1cells)),{'LineWidth',1.6,'Color',[0 0 0],'LineStyle','-'},1);
%hold on;t2=shadedErrorBar(((-Window):1*edgestep:(Window-1))/(samplingf/1000),nanmean(V1selslc),1.96*nanstd(V1selslc)/sqrt(length(V1cells)),{'LineWidth',1.6,'Color',[0 0 1]},1);
xlim([-400 1000]);ylim([-0.5 0.5])
axx=gca;
hold on; patch([0 500 500 0],[axx.YLim(1) axx.YLim(1) axx.YLim(2) axx.YLim(2)],[1 0.7 0.7],...
    'FaceAlpha',0.1,'EdgeColor','none')
xlabel('time (ms)');

 lgd = legend([t1.mainLine t2.mainLine],'baseline selectivity index, correct trials','baseline selectivity index, grooming trials',...
     'Location','southeast');
 lgd.FontSize = 7;



