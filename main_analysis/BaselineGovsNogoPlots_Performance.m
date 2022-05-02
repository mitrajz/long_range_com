
causalbinning = 2; 
                   
frnormalize =0; 
onlyresp = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For this part, I pool FB and FF
Correct = load('OnlyBs_AllBsAnimals_cells_onlyC_40ms_withlaAbs.mat');
Incorrect = load('OnlyBs_AllBsAnimals_cells_onlyIC_40ms_withlaAbs.mat');
Correct = load('cells_onlyC_40ms.mat');
Incorrect = load('cells_onlyIC_40ms.mat');

edgestep = Correct.edgestep;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if onlyresp
    during = floor(size(Correct.V1cells{1}.nbs.go,2)/2):floor(3*size(Correct.V1cells{1}.nbs.go,2)/4);  %200:300
    before =  floor(size(Correct.V1cells{1}.nbs.go,2)/4):floor(size(Correct.V1cells{1}.nbs.go,2)/2);
    a=cellfun(@(x) min(ranksum(sum(x.nbs.go(:,during),2),sum(x.nbs.go(:,before),2)),ranksum(sum(x.nbs.nogo(:,during),2),sum(x.nbs.nogo(:,before),2))),Correct.V1cells);
    in=find(a<0.01);
    b=cellfun(@(x) min((sum(sum(x.nbs.go(:,during)))-sum(sum(x.nbs.go(:,before)))),(sum(sum(x.nbs.nogo(:,during)))-sum(sum(x.nbs.nogo(:,before))))),Correct.V1cells);
    in2=find(b>0);
    Correct.V1cells=Correct.V1cells(intersect(in,in2));
    
    a=cellfun(@(x) min(ranksum(sum(x.nbs.go(:,during),2),sum(x.nbs.go(:,before),2)),ranksum(sum(x.nbs.nogo(:,during),2),sum(x.nbs.nogo(:,before),2))),Correct.LMcells);
    in=find(a<0.01);
    b=cellfun(@(x) min((sum(sum(x.nbs.go(:,during)))-sum(sum(x.nbs.go(:,before)))),(sum(sum(x.nbs.nogo(:,during)))-sum(sum(x.nbs.nogo(:,before))))),Correct.LMcells);
    in2=find(b>0);
    Correct.LMcells=Correct.LMcells(intersect(in,in2));
    
    
    
    during = floor(size(Incorrect.V1cells{1}.nbs.go,2)/2):floor(3*size(Incorrect.V1cells{1}.nbs.go,2)/4);  %200:300
    before =  floor(size(Incorrect.V1cells{1}.nbs.go,2)/4):floor(size(Incorrect.V1cells{1}.nbs.go,2)/2);
    a=cellfun(@(x) min(ranksum(sum(x.nbs.go(:,during),2),sum(x.nbs.go(:,before),2)),ranksum(sum(x.nbs.nogo(:,during),2),sum(x.nbs.nogo(:,before),2))),Incorrect.V1cells);
    in=find(a<0.01);
    b=cellfun(@(x) min((sum(sum(x.nbs.go(:,during)))-sum(sum(x.nbs.go(:,before)))),(sum(sum(x.nbs.nogo(:,during)))-sum(sum(x.nbs.nogo(:,before))))),Incorrect.V1cells);
    in2=find(b>0);
    Incorrect.V1cells=Incorrect.V1cells(intersect(in,in2));
    
    a=cellfun(@(x) min(ranksum(sum(x.nbs.go(:,during),2),sum(x.nbs.go(:,before),2)),ranksum(sum(x.nbs.nogo(:,during),2),sum(x.nbs.nogo(:,before),2))),Incorrect.LMcells);
    in=find(a<0.01);
    b=cellfun(@(x) min((sum(sum(x.nbs.go(:,during)))-sum(sum(x.nbs.go(:,before)))),(sum(sum(x.nbs.nogo(:,during)))-sum(sum(x.nbs.nogo(:,before))))),Incorrect.LMcells);
    in2=find(b>0);
    Incorrect.LMcells=Incorrect.LMcells(intersect(in,in2));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timeoffset = 0;
if frnormalize
    Correct_allV1go=cell2mat(transpose(cellfun(@(x) x.nbsAv.go/x.Avfr ,Correct.V1cells,'UniformOutput',0)));
    Correct_allV1nogo=cell2mat(transpose(cellfun(@(x) x.nbsAv.nogo/x.Avfr ,Correct.V1cells,'UniformOutput',0)));
    Incorrect_allV1go=cell2mat(transpose(cellfun(@(x) x.nbsAv.go/x.Avfr ,Incorrect.V1cells,'UniformOutput',0)));
    Incorrect_allV1nogo=cell2mat(transpose(cellfun(@(x) x.nbsAv.nogo/x.Avfr ,Incorrect.V1cells,'UniformOutput',0)));
else
    Correct_allV1go=cell2mat(transpose(cellfun(@(x) x.nbsAv.go ,Correct.V1cells,'UniformOutput',0)));
    Correct_allV1nogo=cell2mat(transpose(cellfun(@(x) x.nbsAv.nogo ,Correct.V1cells,'UniformOutput',0)));
    Incorrect_allV1go=cell2mat(transpose(cellfun(@(x) x.nbsAv.go ,Incorrect.V1cells,'UniformOutput',0)));
    Incorrect_allV1nogo=cell2mat(transpose(cellfun(@(x) x.nbsAv.nogo ,Incorrect.V1cells,'UniformOutput',0)));
end
if causalbinning == 1
    Correct_allV1go = circshift(Correct_allV1go,1,2);
    Correct_allV1nogo = circshift(Correct_allV1nogo,1,2);
    Incorrect_allV1go = circshift(Incorrect_allV1go,1,2);
    Incorrect_allV1nogo = circshift(Incorrect_allV1nogo,1,2);
elseif causalbinning == 2
    timeoffset = .5*edgestep/(samplingf/1000);
end
figure;plot(timeoffset + ((-Window):1*edgestep:(Window-1))/(samplingf/1000),1000*nanmean(Correct_allV1go,1)/(edgestep/30),'Color',[0 0.8 0],'LineWidth',1.6);
hold on;plot(timeoffset +((-Window):1*edgestep:(Window-1))/(samplingf/1000),1000*nanmean(Correct_allV1nogo,1)/(edgestep/30),'Color',[1 0 0],'LineWidth',1.6);

hold on;plot(timeoffset +((-Window):1*edgestep:(Window-1))/(samplingf/1000),1000*nanmean(Incorrect_allV1go,1)/(edgestep/30),'Color',1*[0 0.8 0],'LineWidth',1.6,'LineStyle','-.');
hold on;plot(timeoffset +((-Window):1*edgestep:(Window-1))/(samplingf/1000),1000*nanmean(Incorrect_allV1nogo,1)/(edgestep/30),'Color',1*[1 0 0],'LineWidth',1.6,'LineStyle','-.');
legend({'Hit','CR','Miss','FA'})
%legend({'Hit','CR','go - grooming','nogo - grooming'})

title('V1 baseline: go (green) vs. nogo (red)')
xlabel('Time (ms)');ylabel('Firing rate (Hz)');
xlim([-100 1000]);ylim([0 18]);
axx = gca;
hold on; patch([0 500 500 0],[axx.YLim(1) axx.YLim(1) axx.YLim(2) axx.YLim(2)],[0.7 0.7 0.7],...
    'FaceAlpha',0.3,'EdgeColor','none')
axx.XTick=0:100:1000;
axx.XTickLabel={'0','','200','','400','','600','','800','','1000'};
axx.YTick=0:5:20;
axx.YTickLabel={'0','5','10','15'};
set(gcf,'Color','w')
legend('boxoff')
%% LM 
timeoffset = 0;
if frnormalize
    Correct_allLMgo=cell2mat(transpose(cellfun(@(x) x.nbsAv.go/x.Avfr ,Correct.LMcells,'UniformOutput',0)));
    Correct_allLMnogo=cell2mat(transpose(cellfun(@(x) x.nbsAv.nogo/x.Avfr ,Correct.LMcells,'UniformOutput',0)));
    Incorrect_allLMgo=cell2mat(transpose(cellfun(@(x) x.nbsAv.go/x.Avfr ,Incorrect.LMcells,'UniformOutput',0)));
    Incorrect_allLMnogo=cell2mat(transpose(cellfun(@(x) x.nbsAv.nogo/x.Avfr ,Incorrect.LMcells,'UniformOutput',0)));
else
    Correct_allLMgo=cell2mat(transpose(cellfun(@(x) x.nbsAv.go ,Correct.LMcells,'UniformOutput',0)));
    Correct_allLMnogo=cell2mat(transpose(cellfun(@(x) x.nbsAv.nogo ,Correct.LMcells,'UniformOutput',0)));
    Incorrect_allLMgo=cell2mat(transpose(cellfun(@(x) x.nbsAv.go ,Incorrect.LMcells,'UniformOutput',0)));
    Incorrect_allLMnogo=cell2mat(transpose(cellfun(@(x) x.nbsAv.nogo ,Incorrect.LMcells,'UniformOutput',0)));
end
if causalbinning == 1
    Correct_allLMgo = circshift(Correct_allLMgo,1,2);
    Correct_allLMnogo = circshift(Correct_allLMnogo,1,2);
    Incorrect_allLMgo = circshift(Incorrect_allLMgo,1,2);
    Incorrect_allLMnogo = circshift(Incorrect_allLMnogo,1,2);
elseif causalbinning == 2
    timeoffset = .5*edgestep/(samplingf/1000);
end
figure;plot(timeoffset +((-Window):1*edgestep:(Window-1))/(samplingf/1000),1000*nanmean(Correct_allLMgo)/(edgestep/30),'Color',[0 0.8 0],'LineWidth',1.6);
hold on;plot(timeoffset +((-Window):1*edgestep:(Window-1))/(samplingf/1000),1000*nanmean(Correct_allLMnogo)/(edgestep/30),'Color',[1 0 0],'LineWidth',1.6);
hold on;plot(timeoffset +((-Window):1*edgestep:(Window-1))/(samplingf/1000),1000*nanmean(Incorrect_allLMgo)/(edgestep/30),'Color',1*[0 0.8 0],'LineWidth',1.6,'LineStyle','-.');
hold on;plot(timeoffset +((-Window):1*edgestep:(Window-1))/(samplingf/1000),1000*nanmean(Incorrect_allLMnogo)/(edgestep/30),'Color',1*[1 0 0],'LineWidth',1.6,'LineStyle','-.');
legend({'Hit','CR','Miss','FA'})

title('LM baseline: go (green) vs. nogo (red)')
xlabel('Time (ms)');ylabel('Firing rate (Hz)');
xlim([-100 1000]);ylim([0 18]);
axx = gca;
hold on; patch([0 500 500 0],[axx.YLim(1) axx.YLim(1) axx.YLim(2) axx.YLim(2)],[0.7 0.7 0.7],...
    'FaceAlpha',0.3,'EdgeColor','none')
axx.XTick=0:100:1000;
axx.XTickLabel={'0','','200','','400','','600','','800','','1000'};
axx.YTick=0:5:20;
axx.YTickLabel={'0','5','10','15'};
set(gcf,'Color','w')
legend('boxoff')