function BaselineGovsNogoPlots(V1Cells,LMCells,Window,edgestep,samplingf,frnormalize,causalbinning,BehCells)
stimlengthinms = 500; %500
zoomedin=0;
dodif= 0;
doerbars = 1
%% raw traces og LM and V1 (all units) in go vs. nogo
% For this part, I pool FB and FF
% 10 ms bins might be better
if frnormalize
    allV1go=cell2mat(transpose(cellfun(@(x) x.nbsAv.go/x.Avfr ,V1Cells,'UniformOutput',0)));
    allV1nogo=cell2mat(transpose(cellfun(@(x) x.nbsAv.nogo/x.Avfr ,V1Cells,'UniformOutput',0)));
else
%     allV1go=cell2mat(transpose(cellfun(@(x) x.nbsAv.go ,V1Cells(find(cellfun(@(x) size(x.nbs.go,1),V1Cells)>5)),'UniformOutput',0)));
%     allV1nogo=cell2mat(transpose(cellfun(@(x) x.nbsAv.nogo ,V1Cells(find(cellfun(@(x) size(x.nbs.nogo,1),V1Cells)>5)),'UniformOutput',0)));
     allV1go=cell2mat(transpose(cellfun(@(x) x.nbsAv.go ,V1Cells,'UniformOutput',0)));
     allV1nogo=cell2mat(transpose(cellfun(@(x) x.nbsAv.nogo ,V1Cells,'UniformOutput',0)));

end
timeoffset = 0;
if causalbinning == 1
    allV1go = circshift(allV1go,1,2);
    allV1nogo = circshift(allV1nogo,1,2);
elseif causalbinning == 2
    timeoffset = .5*edgestep/(samplingf/1000);
end
figure;plot(timeoffset + ((-Window):1*edgestep:(Window-1))/(samplingf/1000),1000*nanmean(allV1go,1)/(edgestep/30),'Color',[0 0.8 0],'LineWidth',1.6);
hold on;plot(timeoffset + ((-Window):1*edgestep:(Window-1))/(samplingf/1000),1000*nanmean(allV1nogo,1)/(edgestep/30),'Color',[1 0 0],'LineWidth',1.6);

if doerbars
    allV1go_rand_mean=nan(1000,size(allV1go,2));
    allV1nogo_rand_mean=nan(1000,size(allV1nogo,2));
    for randsamp = 1:1000
        allV1go_rand = cell2mat(transpose(cellfun(@(x) nanmean(x.nbs.go(randi(size(x.nbs.go,1),1,size(x.nbs.go,1)),:)),...
            V1Cells,'UniformOutput',0)));
        allV1nogo_rand = cell2mat(transpose(cellfun(@(x) nanmean(x.nbs.nogo(randi(size(x.nbs.nogo,1),1,size(x.nbs.nogo,1)),:)),...
            V1Cells,'UniformOutput',0)));
        allV1go_rand_mean(randsamp,:) = nanmean(allV1go_rand,1)/(edgestep/30);
        allV1nogo_rand_mean(randsamp,:) = nanmean(allV1nogo_rand,1)/(edgestep/30);
    end
    hold on;patch([timeoffset + ((-Window):1*edgestep:(Window-1))/(samplingf/1000),fliplr(timeoffset + ((-Window):1*edgestep:(Window-1))/(samplingf/1000))],...
    [1000*prctile(allV1go_rand_mean,97.5),fliplr(1000*prctile(allV1go_rand_mean,2.5))],[0 0.8 0],'FaceAlpha',0.2,'EdgeAlpha',0);
    hold on;patch([timeoffset + ((-Window):1*edgestep:(Window-1))/(samplingf/1000),fliplr(timeoffset + ((-Window):1*edgestep:(Window-1))/(samplingf/1000))],...
    [1000*prctile(allV1nogo_rand_mean,97.5),fliplr(1000*prctile(allV1nogo_rand_mean,2.5))],[1 0 0],'FaceAlpha',0.2,'EdgeAlpha',0);


end


%title('V1 baseline: go (green) vs. nogo (red)')
xlabel('Time (ms)');ylabel('Firing rate (Hz)');
xlim([-100 1000]);ylim([0 18]);
if zoomedin
    ylim([0 25]);
end
axx = gca;
hold on; patch([0 stimlengthinms stimlengthinms 0],[axx.YLim(1) axx.YLim(1) axx.YLim(2) axx.YLim(2)],[0.7 0.7 0.7],...
    'FaceAlpha',0.3,'EdgeColor','none')
axx.XTick=0:100:1000;
axx.XTickLabel={'0','','200','','400','','600','','800','','1000'};
if zoomedin
    xlim([0 200]);
    axx.XTick=0:50:200;
    axx.XTickLabel={'0','','100','','200'};
end
axx.YTick=0:5:20;
axx.YTickLabel={'0','5','10','15','20','25'};
set(gcf,'Color','w')
legend({'go - correct','no go - correct'})
legend('boxoff')
box(axx,'off');
if zoomedin
    legend('hide');
end

if frnormalize
    allLMgo=cell2mat(transpose(cellfun(@(x) x.nbsAv.go/x.Avfr ,LMCells,'UniformOutput',0)));
    allLMnogo=cell2mat(transpose(cellfun(@(x) x.nbsAv.nogo/x.Avfr ,LMCells,'UniformOutput',0)));
else
    allLMgo=cell2mat(transpose(cellfun(@(x) x.nbsAv.go ,LMCells,'UniformOutput',0)));
    allLMnogo=cell2mat(transpose(cellfun(@(x) x.nbsAv.nogo ,LMCells,'UniformOutput',0)));
end
timeoffset = 0;
if causalbinning == 1
    allLMgo = circshift(allLMgo,1,2);
    allLMnogo = circshift(allLMnogo,1,2);
elseif causalbinning == 2
    timeoffset = .5*edgestep/(samplingf/1000);
end
figure;plot(timeoffset + ((-Window):1*edgestep:(Window-1))/(samplingf/1000),1000*nanmean(allLMgo)/(edgestep/30),'Color',[0 0.8 0],'LineWidth',1.6);
hold on;plot(timeoffset + ((-Window):1*edgestep:(Window-1))/(samplingf/1000),1000*nanmean(allLMnogo)/(edgestep/30),'Color',[1 0 0],'LineWidth',1.6);


%title('LM baseline: go (green) vs. nogo (red)')
xlabel('Time (ms)');ylabel('Firing rate (Hz)');
xlim([-100 1000]);ylim([0 18]);
if zoomedin
    ylim([0 25]);
end
axx = gca;
hold on; patch([0 stimlengthinms stimlengthinms 0],[axx.YLim(1) axx.YLim(1) axx.YLim(2) axx.YLim(2)],[0.7 0.7 0.7],...
    'FaceAlpha',0.3,'EdgeColor','none')
axx.XTick=0:100:1000;
axx.XTickLabel={'0','','200','','400','','600','','800','','1000'};
if zoomedin
    xlim([0 200]);
    axx.XTick=0:50:200;
    axx.XTickLabel={'0','','100','','200'};
end
axx.YTick=0:5:20;
axx.YTickLabel={'0','5','10','15','20','25'};
set(gcf,'Color','w')
legend({'go - correct','no go - correct'})
legend('boxoff')
box(axx,'off');
if zoomedin
    legend('hide');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% licks
alllicksgo=cell2mat(cellfun(@(x) nanmean(x.licks.go,1),BehCells,'UniformOutput',0)');
alllicksnogo=cell2mat(cellfun(@(x) nanmean(x.licks.nogo,1),BehCells,'UniformOutput',0)');

timeoffset = 0;
if causalbinning == 1
    alllicksgo = circshift(alllicksgo,1,2);
    alllicksnogo = circshift(alllicksnogo,1,2);
elseif causalbinning == 2
    timeoffset = .5*edgestep/(samplingf/1000);
end
figure;plot(timeoffset + ((-Window):1*edgestep:(Window-1))/(samplingf/1000),nanmean(alllicksgo,1)/(edgestep),'Color',[0 0.8 0],'LineWidth',1.6);
hold on;plot(timeoffset + ((-Window):1*edgestep:(Window-1))/(samplingf/1000),nanmean(alllicksnogo,1)/(edgestep),'Color',[1 0 0],'LineWidth',1.6);
%title('licks: go (green) vs. nogo (red)')
xlabel('Time (ms)');
xlim([-100 1000]);%ylim([0 25]);
axx = gca;
hold on; patch([0 stimlengthinms stimlengthinms 0],[axx.YLim(1) axx.YLim(1) axx.YLim(2) axx.YLim(2)],[0.7 0.7 0.7],...
    'FaceAlpha',0.3,'EdgeColor','none')
axx.XTick=0:100:1000;
set(gcf,'Color','w')
legend({'go - correct','no go - correct'})
legend('boxoff')

%% imagesc : 2 figures: sorted off go or nogo responsivity
clm=0.5;

[~,sorti_V1_go] = sort(sum(allV1go(:,floor(size(allV1go,2)/2):end),2) - sum(allV1go(:,1:floor(size(allV1go,2)/2)),2));
[~,sorti_V1_nogo] = sort(sum(allV1nogo(:,floor(size(allV1nogo,2)/2):end),2) - sum(allV1nogo(:,1:floor(size(allV1nogo,2)/2)),2));
[~,sorti_LM_go] = sort(sum(allLMgo(:,floor(size(allLMgo,2)/2):end),2) - sum(allLMgo(:,1:floor(size(allLMgo,2)/2)),2));
[~,sorti_LM_nogo] = sort(sum(allLMnogo(:,floor(size(allLMnogo,2)/2):end),2) - sum(allLMnogo(:,1:floor(size(allLMnogo,2)/2)),2));
%%%%%%% Cells sorted based on go response size
indv_sorting = 1;
sortingtype = 'go'; % only counts when indv_sorting = 0


f = figure;
f.Units = 'normalized';f.Position = [0.25 0.25 0.5 0.4];
if ~indv_sorting
    f.Name = sprintf('V1 Cells sorted based on %s responses',sortingtype);
    s1=subplot(1,2,1);imagesc(timeoffset + ((-Window):1*edgestep:(Window-1))/(samplingf/1000),[],allV1go(eval(sprintf('sorti_V1_%s',sortingtype)),:));
else
    f.Name = 'V1 Cells sorted separately for different stimuli';
    s1=subplot(1,2,1);imagesc(timeoffset + ((-Window):1*edgestep:(Window-1))/(samplingf/1000),[],allV1go(eval(sprintf('sorti_V1_%s','go')),:));
end
greymap=colormap('gray');
colormap(flipud(greymap));
set(gca,'CLim',[0 clm]);
xlim([-100 1000]);
s1.Title.String = 'V1 - go';
s1.Title.FontWeight = 'normal';
s1.Title.FontSize = 9;
s1.FontSize=8;V1Cells
s1.XLabel.String = 'Time (ms)';
s1.YLabel.String = 'Cells';


if ~indv_sorting
    s2=subplot(1,2,2);imagesc(timeoffset + ((-Window):1*edgestep:(Window-1))/(samplingf/1000),[],allV1nogo(eval(sprintf('sorti_V1_%s',sortingtype)),:));
else
    s2=subplot(1,2,2);imagesc(timeoffset + ((-Window):1*edgestep:(Window-1))/(samplingf/1000),[],allV1nogo(eval(sprintf('sorti_V1_%s','nogo')),:));
end
greymap=colormap('gray');
colormap(flipud(greymap));
set(gca,'CLim',[0 clm]);
xlim([-100 1000]);
s2.Title.String = 'V1 - no go - correct';
s2.Title.FontWeight = 'normal';
s2.Title.FontSize = 9;
s2.FontSize=8;
s2.XLabel.String = 'Time (ms)';
s2.YLabel.String = 'Cells';

f = figure;
f.Units = 'normalized';f.Position = [0.25 0.25 0.5 0.4];
if ~indv_sorting
    f.Name = sprintf('V1 Cells sorted based on %s responses',sortingtype);
    s3=subplot(1,2,1);imagesc(timeoffset + ((-Window):1*edgestep:(Window-1))/(samplingf/1000),[],allLMgo(eval(sprintf('sorti_LM_%s',sortingtype)),:));
else
    f.Name = 'V1 Cells sorted separately for different stimuli';
    s3=subplot(1,2,1);imagesc(timeoffset + ((-Window):1*edgestep:(Window-1))/(samplingf/1000),[],allLMgo(eval(sprintf('sorti_LM_%s','go')),:));
end
greymap=colormap('gray');
colormap(flipud(greymap));
set(gca,'CLim',[0 clm]);
xlim([-100 1000]);
s3.Title.String = 'LM - go';
s3.Title.FontWeight = 'normal';
s3.Title.FontSize = 9;
s3.FontSize=8;
s3.XLabel.String = 'Time (ms)';
s3.YLabel.String = 'Cells';

if ~indv_sorting
    s4=subplot(1,2,2);imagesc(timeoffset + ((-Window):1*edgestep:(Window-1))/(samplingf/1000),[],allLMnogo(eval(sprintf('sorti_LM_%s',sortingtype)),:));
else
    s4=subplot(1,2,2);imagesc(timeoffset + ((-Window):1*edgestep:(Window-1))/(samplingf/1000),[],allLMnogo(eval(sprintf('sorti_LM_%s','nogo')),:));
end
greymap=colormap('gray');
colormap(flipud(greymap));
set(gca,'CLim',[0 clm]);
xlim([-100 1000]);
s4.Title.String = 'LM - no go - correct';
s4.Title.FontWeight = 'normal';
s4.Title.FontSize = 9;
s4.FontSize=8;
s4.XLabel.String = 'Time (ms)';
s4.YLabel.String = 'Cells';

set(gcf,'Color','w');
%% difference plots - Version1
if dodif
V1_pool = [allV1nogo;allV1go];
LM_pool = [allLMnogo;allLMgo];

LMrand=nan(50000,size(allLMgo,2));
V1rand=nan(50000,size(allV1go,2));

for rrand = 1:50000
    selectV1=zeros(1,size(allV1go,1)+size(allV1nogo,1))';
    selectV1(randi(length(selectV1),[1,size(allV1nogo,1)]))=1;
    V1rand(rrand,:) = 1000*(nanmean(V1_pool((find(~selectV1)),:))-nanmean(V1_pool((find(selectV1)),:)))/(edgestep/30);
    
    selectLM=zeros(1,size(allLMgo,1)+size(allLMnogo,1))';
    selectLM(randi(length(selectLM),[1,size(allLMnogo,1)]))=1;
    LMrand(rrand,:) = 1000*(nanmean(LM_pool((find(~selectLM)),:))-nanmean(LM_pool((find(selectLM)),:)))/(edgestep/30);
end



figure;plot(timeoffset + ((-Window):1*edgestep:(Window-1))/(samplingf/1000),1000*(nanmean(allV1go,1)-nanmean(allV1nogo,1))/(edgestep/30),'Color','b','LineWidth',1.6);
xpoints=[timeoffset + ((-Window):1*edgestep:(Window-1))/(samplingf/1000) fliplr(timeoffset + ((-Window):1*edgestep:(Window-1))/(samplingf/1000))];
ypoints=[prctile(V1rand,5),fliplr(prctile(V1rand,95))];
hold on;patch(xpoints,ypoints,'b','FaceAlpha',0.2,'EdgeAlpha',0)

%title('V1 baseline: go - nogo')
xlabel('Time (ms)');ylabel('Firing rate (Hz)');
xlim([-100 1000]);ylim([-20 20]);
axx = gca;
hold on; patch([0 stimlengthinms stimlengthinms 0],[axx.YLim(1) axx.YLim(1) axx.YLim(2) axx.YLim(2)],[0.7 0.7 0.7],...
    'FaceAlpha',0.3,'EdgeColor','none')
set(gcf,'Color','w')
axx.XTick=0:100:1000;
axx.XTickLabel={'0','','200','','400','','600','','800','','1000'};
legend({'V1 (go - nogo)'})
legend('boxoff')


figure;plot(timeoffset + ((-Window):1*edgestep:(Window-1))/(samplingf/1000),1000*(nanmean(allLMgo)-nanmean(allLMnogo))/(edgestep/30),'Color','b','LineWidth',1.6);
xpoints=[timeoffset + ((-Window):1*edgestep:(Window-1))/(samplingf/1000) fliplr(timeoffset + ((-Window):1*edgestep:(Window-1))/(samplingf/1000))];
ypoints=[prctile(LMrand,5),fliplr(prctile(LMrand,95))];
hold on;patch(xpoints,ypoints,'b','FaceAlpha',0.2,'EdgeAlpha',0)

%title('LM baseline: go - nogo')
xlabel('Time (ms)');ylabel('Firing rate (Hz)');
xlim([-100 1000]);ylim([-10 10]);
axx = gca;
hold on; patch([0 stimlengthinms stimlengthinms 0],[axx.YLim(1) axx.YLim(1) axx.YLim(2) axx.YLim(2)],[0.7 0.7 0.7],...
    'FaceAlpha',0.3,'EdgeColor','none')
set(gcf,'Color','w')
axx.XTick=0:100:1000;
axx.XTickLabel={'0','','200','','400','','600','','800','','1000'};
legend({'LM (go - nogo)'})
legend('boxoff')
%% difference plots - Version2
selindex=1;

LMrand=nan(1000,size(LMCells{1}.nbs.go,2));% 1000* time
V1rand=nan(1000,size(V1Cells{1}.nbs.go,2));% 1000 * time
% for V1
for rrand = 1:1000  
    onerep_gonodif_percell_overtime_V1 = nan(numel(V1Cells),size(V1Cells{1}.nbs.go,2));
    for i= 1:numel(V1Cells)
        
        V1_pool = [V1Cells{i}.nbs.go;V1Cells{i}.nbs.nogo];
        % LM_pool = [LMCells{i}.nbs.go;LMCells{i}.nbs.nogo];
        
        selectV1=zeros(1,size(V1Cells{i}.nbs.go,1)+size(V1Cells{i}.nbs.nogo,1))';
        selectV1(randi(length(selectV1),[1,size(V1Cells{i}.nbs.nogo,1)]))=1;
        % gonogodif percell
        if selindex
            onerep_gonodif_percell_overtime_V1(i,:) = ((nanmean(V1_pool((find(~selectV1)),:))-nanmean(V1_pool((find(selectV1)),:)))./ ...
            (nanmean(V1_pool((find(~selectV1)),:))+nanmean(V1_pool((find(selectV1)),:))) )...
            ;
        else
            onerep_gonodif_percell_overtime_V1(i,:) = (nanmean(V1_pool((find(~selectV1)),:))-nanmean(V1_pool((find(selectV1)),:)))...
            ;
        end
    end
    V1rand(rrand,:) = nanmean(zscore(onerep_gonodif_percell_overtime_V1')',1) ;
end
% for LM
for rrand = 1:1000   
    onerep_gonodif_percell_overtime_LM = nan(numel(LMCells),size(LMCells{1}.nbs.go,2));
    for i= 1:numel(LMCells)
        
        LM_pool = [LMCells{i}.nbs.go;LMCells{i}.nbs.nogo];
        % LM_pool = [LMCells{i}.nbs.go;LMCells{i}.nbs.nogo];
        
        selectLM=zeros(1,size(LMCells{i}.nbs.go,1)+size(LMCells{i}.nbs.nogo,1))';
        selectLM(randi(length(selectLM),[1,size(LMCells{i}.nbs.nogo,1)]))=1;
        % gonogodif percell
        if selindex
            onerep_gonodif_percell_overtime_LM(i,:) = ((nanmean(LM_pool((find(~selectLM)),:))-nanmean(LM_pool((find(selectLM)),:)))./ ...
            (nanmean(LM_pool((find(~selectLM)),:))+nanmean(LM_pool((find(selectLM)),:))) )...
            ;
        else
            onerep_gonodif_percell_overtime_LM(i,:) = (nanmean(LM_pool((find(~selectLM)),:))-nanmean(LM_pool((find(selectLM)),:)))...
            ;
        end
    end
    LMrand(rrand,:) = nanmean(zscore(onerep_gonodif_percell_overtime_LM')',1);
end
%%
if selindex
    gonogodifpercell_V1 = (allV1go-allV1nogo)./(allV1go+allV1nogo); 
else
    gonogodifpercell_V1 = (allV1go-allV1nogo); 
end
figure;plot(timeoffset + ((-Window):1*edgestep:(Window-1))/(samplingf/1000),(nanmean(zscore(gonogodifpercell_V1')',1)),'Color','b','LineWidth',1.6);
xpoints=[timeoffset + ((-Window):1*edgestep:(Window-1))/(samplingf/1000) fliplr(timeoffset + ((-Window):1*edgestep:(Window-1))/(samplingf/1000))];
ypoints=[prctile(V1rand,2.5),fliplr(prctile(V1rand,97.5))];
hold on;patch(xpoints,ypoints,'b','FaceAlpha',0.2,'EdgeAlpha',0)

%title('V1 baseline: go - nogo')
xlabel('Time (ms)');
%ylabel('Firing rate (Hz)');
xlim([-100 1000]);ylim([-1 1]);
axx = gca;
hold on; patch([0 stimlengthinms stimlengthinms 0],[axx.YLim(1) axx.YLim(1) axx.YLim(2) axx.YLim(2)],[0.7 0.7 0.7],...
    'FaceAlpha',0.3,'EdgeColor','none')
set(gcf,'Color','w')
axx.XTick=0:100:1000;
axx.XTickLabel={'0','','200','','400','','600','','800','','1000'};
% legend({'V1 (go - nogo)'})
% legend('boxoff')


if selindex
    gonogodifpercell_LM = (allLMgo-allLMnogo)./(allLMgo+allLMnogo); 
else
    gonogodifpercell_LM = (allLMgo-allLMnogo); 
end

figure;plot(timeoffset + ((-Window):1*edgestep:(Window-1))/(samplingf/1000),(nanmean(zscore(gonogodifpercell_LM')',1)),'Color','b','LineWidth',1.6);
xpoints=[timeoffset + ((-Window):1*edgestep:(Window-1))/(samplingf/1000) fliplr(timeoffset + ((-Window):1*edgestep:(Window-1))/(samplingf/1000))];
ypoints=[prctile(LMrand,2.5),fliplr(prctile(LMrand,97.5))];
hold on;patch(xpoints,ypoints,'b','FaceAlpha',0.2,'EdgeAlpha',0)

%title('LM baseline: go - nogo')
xlabel('Time (ms)');ylabel('Firing rate (Hz)');
xlim([-100 1000]);ylim([-1 1]);
axx = gca;
hold on; patch([0 stimlengthinms stimlengthinms 0],[axx.YLim(1) axx.YLim(1) axx.YLim(2) axx.YLim(2)],[0.7 0.7 0.7],...
    'FaceAlpha',0.3,'EdgeColor','none')
set(gcf,'Color','w')
axx.XTick=0:100:1000;
axx.XTickLabel={'0','','200','','400','','600','','800','','1000'};
legend({'LM (go - nogo)'})
legend('boxoff')
end