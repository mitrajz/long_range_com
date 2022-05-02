
frnormalize = 0 ;
causalbinning = 2; 
V1selectsign = 0; % if 0: alltogether, 1: positive, 2:negative
LMselectsign = 0; % if 0: alltogether, 1: positive, 2:negative

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% V1
figure;
for lag=1:8
subplot(8,1,lag)
timeoffset = 0;

V1signgoperlag = cellfun(@(x) x.smb.go(lag) < 0, V1cells);
V1signnogoperlag = cellfun(@(x) x.smb.nogo(lag) < 0, V1cells);


if frnormalize
    allV1go=cell2mat(transpose(cellfun(@(x) x.nlsAv.go(lag,:)/x.Avfr ,V1cells,'UniformOutput',0)));
    allV1nogo=cell2mat(transpose(cellfun(@(x) x.nlsAv.nogo(lag,:)/x.Avfr ,V1cells,'UniformOutput',0)));
else
    allV1go=cell2mat(transpose(cellfun(@(x) x.nlsAv.go(lag,:) ,V1cells,'UniformOutput',0)));
    allV1nogo=cell2mat(transpose(cellfun(@(x) x.nlsAv.nogo(lag,:) ,V1cells,'UniformOutput',0)));
end
if causalbinning == 1
    allV1go = circshift(allV1go,1,2);
    allV1nogo = circshift(allV1nogo,1,2);
elseif causalbinning == 2
    timeoffset = .5*edgestep/(samplingf/1000);
end

if V1selectsign == 1
    allV1go = allV1go(find(V1signgoperlag == 1),:);
    allV1nogo = allV1nogo(find(V1signnogoperlag == 1),:);
elseif V1selectsign == 2
    allV1go = allV1go(find(V1signgoperlag == 0),:);
    allV1nogo = allV1nogo(find(V1signnogoperlag == 0),:);
end

plot(timeoffset+((-Window):1*edgestep:(Window-1))/(samplingf/1000),1000*nanmean(allV1go,1)/(edgestep/30),'Color',[0 0.8 0],'LineWidth',1.6,'LineStyle','--');
hold on;plot(timeoffset+((-Window):1*edgestep:(Window-1))/(samplingf/1000),1000*nanmean(allV1nogo,1)/(edgestep/30),'Color',[1 0 0],'LineWidth',1.6,'LineStyle','--');



hold on;
timeoffset = 0;

if frnormalize
    allV1go=cell2mat(transpose(cellfun(@(x) x.nbsAv.go/x.Avfr ,V1cells,'UniformOutput',0)));
    allV1nogo=cell2mat(transpose(cellfun(@(x) x.nbsAv.nogo/x.Avfr ,V1cells,'UniformOutput',0)));
else
    allV1go=cell2mat(transpose(cellfun(@(x) x.nbsAv.go ,V1cells,'UniformOutput',0)));
    allV1nogo=cell2mat(transpose(cellfun(@(x) x.nbsAv.nogo ,V1cells,'UniformOutput',0)));
end
if causalbinning == 1
    allV1go = circshift(allV1go,1,2);
    allV1nogo = circshift(allV1nogo,1,2);
elseif causalbinning == 2
    timeoffset = .5*edgestep/(samplingf/1000);
end


if V1selectsign == 1
    allV1go = allV1go(find(V1signgoperlag == 1),:);
    allV1nogo = allV1nogo(find(V1signnogoperlag == 1),:);
elseif V1selectsign == 2
    allV1go = allV1go(find(V1signgoperlag == 0),:);
    allV1nogo = allV1nogo(find(V1signnogoperlag == 0),:);
end


plot(timeoffset+((-Window):1*edgestep:(Window-1))/(samplingf/1000),1000*nanmean(allV1go,1)/(edgestep/30),'Color',[0 0.8 0],'LineWidth',1.6);
hold on;plot(timeoffset+((-Window):1*edgestep:(Window-1))/(samplingf/1000),1000*nanmean(allV1nogo,1)/(edgestep/30),'Color',[1 0 0],'LineWidth',1.6);
if lag == 1
title('V1 baseline: go (green) vs. nogo (red)')
xlabel('time (ms)');
end
if lag == 4
ylabel('average firing rate (Hz)');
end
xlim([-100 1000]);ylim([0 18]);%18
axx = gca;
hold on; patch([0 500 500 0],[axx.YLim(1) axx.YLim(1) axx.YLim(2) axx.YLim(2)],[0.7 0.7 0.7],...
    'FaceAlpha',0.3,'EdgeColor','none')
axx.XTick=0:100:1000;
set(gcf,'Color','w')


startpoint=V1cells{1}.smb_centers.go(lag);
dur=40;
dur2=75;
dur3=150;
% hold on;patch([ startpoint startpoint+dur startpoint+dur startpoint],[0 0 10 10],'b','FaceAlpha',0.1)
hold on;patch([ startpoint startpoint+dur3 startpoint+dur3 startpoint],[0 0 100 100],'b','FaceAlpha',0.1)
hold on;patch([ startpoint startpoint+dur2 startpoint+dur2 startpoint],[0 0 100 100],'b','FaceAlpha',0.1,'EdgeAlpha',0)

set(axx,'Position',[0.13,(0.95-0.1*lag),0.775,0.1])
if lag~=8
    axx.XTickLabel=[];
    axx.YTickLabel=[];
end
end
 %% LM
figure;
for lag=1:8
subplot(8,1,lag)



LMsigngoperlag = cellfun(@(x) x.smb.go(lag) < 0, LMcells);
LMsignnogoperlag = cellfun(@(x) x.smb.nogo(lag) < 0, LMcells);

timeoffset = 0;
if frnormalize
    allLMgo=cell2mat(transpose(cellfun(@(x) x.nlsAv.go(lag,:)/x.Avfr ,LMcells,'UniformOutput',0)));
    allLMnogo=cell2mat(transpose(cellfun(@(x) x.nlsAv.nogo(lag,:)/x.Avfr ,LMcells,'UniformOutput',0)));
else
    allLMgo=cell2mat(transpose(cellfun(@(x) x.nlsAv.go(lag,:) ,LMcells,'UniformOutput',0)));
    allLMnogo=cell2mat(transpose(cellfun(@(x) x.nlsAv.nogo(lag,:) ,LMcells,'UniformOutput',0)));
end
if causalbinning == 1
    allLMgo = circshift(allLMgo,1,2);
    allLMnogo = circshift(allLMnogo,1,2);
elseif causalbinning == 2
    timeoffset = .5*edgestep/(samplingf/1000);

end


if LMselectsign == 1
    allLMgo = allLMgo(find(LMsigngoperlag == 1),:);
    allLMnogo = allLMnogo(find(LMsignnogoperlag == 1),:);
elseif LMselectsign == 2
    allLMgo = allLMgo(find(LMsigngoperlag == 0),:);
    allLMnogo = allLMnogo(find(LMsignnogoperlag == 0),:);
end


plot(timeoffset+((-Window):1*edgestep:(Window-1))/(samplingf/1000),1000*mean(allLMgo)/(edgestep/30),'Color',[0 0.8 0],'LineWidth',1.6,'LineStyle','--');
hold on;plot(timeoffset+((-Window):1*edgestep:(Window-1))/(samplingf/1000),1000*mean(allLMnogo)/(edgestep/30),'Color',[1 0 0],'LineWidth',1.6,'LineStyle','--');



hold on;
timeoffset = 0;
if frnormalize
    allLMgo=cell2mat(transpose(cellfun(@(x) x.nbsAv.go/x.Avfr ,LMcells,'UniformOutput',0)));
    allLMnogo=cell2mat(transpose(cellfun(@(x) x.nbsAv.nogo/x.Avfr ,LMcells,'UniformOutput',0)));
else
    allLMgo=cell2mat(transpose(cellfun(@(x) x.nbsAv.go ,LMcells,'UniformOutput',0)));
    allLMnogo=cell2mat(transpose(cellfun(@(x) x.nbsAv.nogo ,LMcells,'UniformOutput',0)));
end
if causalbinning == 1
    allLMgo = circshift(allLMgo,1,2);
    allLMnogo = circshift(allLMnogo,1,2);
elseif causalbinning == 2
    timeoffset = .5*edgestep/(samplingf/1000);
end

if LMselectsign == 1
    allLMgo = allLMgo(find(LMsigngoperlag == 1),:);
    allLMnogo = allLMnogo(find(LMsignnogoperlag == 1),:);
elseif LMselectsign == 2
    allLMgo = allLMgo(find(LMsigngoperlag == 0),:);
    allLMnogo = allLMnogo(find(LMsignnogoperlag == 0),:);
end


plot(timeoffset+((-Window):1*edgestep:(Window-1))/(samplingf/1000),1000*mean(allLMgo)/(edgestep/30),'Color',[0 0.8 0],'LineWidth',1.6);
hold on;plot(timeoffset+((-Window):1*edgestep:(Window-1))/(samplingf/1000),1000*mean(allLMnogo)/(edgestep/30),'Color',[1 0 0],'LineWidth',1.6);
if lag == 1
title('LM baseline: go (green) vs. nogo (red)')
xlabel('time (ms)');
end
if lag == 4
ylabel('average firing rate (Hz)');
end
xlim([-100 1000]);ylim([0 33]);
axx = gca;
hold on; patch([0 500 500 0],[axx.YLim(1) axx.YLim(1) axx.YLim(2) axx.YLim(2)],[0.7 0.7 0.7],...
    'FaceAlpha',0.3,'EdgeColor','none')
axx.XTick=0:100:1000;
set(gcf,'Color','w')


startpoint=LMcells{1}.smb_centers.go(lag);
% dur=40;
%  hold on;patch([ startpoint startpoint+dur startpoint+dur startpoint],[0 0 10 10],'b','FaceAlpha',0.2)
dur=40;
dur2=75;

dur3=150;
% hold on;patch([ startpoint startpoint+dur startpoint+dur startpoint],[0 0 10 10],'b','FaceAlpha',0.1)
hold on;patch([ startpoint startpoint+dur3 startpoint+dur3 startpoint],[0 0 200 200],'b','FaceAlpha',0.1)
hold on;patch([ startpoint startpoint+dur2 startpoint+dur2 startpoint],[0 0 200 200],'b','FaceAlpha',0.1,'EdgeAlpha',0)

set(axx,'Position',[0.13,(0.95-0.1*lag),0.775,0.1])
if lag~=8
    axx.XTickLabel=[];
    axx.YTickLabel=[];
end

end