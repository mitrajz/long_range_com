function SilencedGovsNogoPlots(V1cells,LMcells,Window,edgestep,samplingf,frnormalize,lag,causalbinning)
%% raw traces og LM and V1 (all units) in go vs. nogo
% For this part, I pool FB and FF
% 10 ms bins might be better
stimlengthinms = 500;
timeoffset=0;
doerbars = 0;
%%
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
figure;plot(timeoffset+((-Window):1*edgestep:(Window-1))/(samplingf/1000),1000*mean(allV1go)/(edgestep/30),'Color',[0 0.8 0],'LineWidth',1.6);
hold on;plot(timeoffset+((-Window):1*edgestep:(Window-1))/(samplingf/1000),1000*mean(allV1nogo/(edgestep/30)),'Color',[1 0 0],'LineWidth',1.6);

if doerbars
    allV1go_rand_mean=nan(1000,size(allV1go,2));
    allV1nogo_rand_mean=nan(1000,size(allV1nogo,2));
    for randsamp = 1:1000
        allV1go_rand = cell2mat(transpose(cellfun(@(x) nanmean(x.nls.go{lag}(randi(size(x.nls.go{lag},1),1,size(x.nls.go{lag},1)),:)),...
            V1cells,'UniformOutput',0)));
        allV1nogo_rand = cell2mat(transpose(cellfun(@(x) nanmean(x.nls.nogo{lag}(randi(size(x.nls.nogo{lag},1),1,size(x.nls.nogo{lag},1)),:)),...
            V1cells,'UniformOutput',0)));
        allV1go_rand_mean(randsamp,:) = nanmean(allV1go_rand,1)/(edgestep/30);
        allV1nogo_rand_mean(randsamp,:) = nanmean(allV1nogo_rand,1)/(edgestep/30);
    end
    hold on;patch([timeoffset + ((-Window):1*edgestep:(Window-1))/(samplingf/1000),fliplr(timeoffset + ((-Window):1*edgestep:(Window-1))/(samplingf/1000))],...
    [1000*prctile(allV1go_rand_mean,97.5),fliplr(1000*prctile(allV1go_rand_mean,2.5))],[0 0.8 0],'FaceAlpha',0.2,'EdgeAlpha',0);
    hold on;patch([timeoffset + ((-Window):1*edgestep:(Window-1))/(samplingf/1000),fliplr(timeoffset + ((-Window):1*edgestep:(Window-1))/(samplingf/1000))],...
    [1000*prctile(allV1nogo_rand_mean,97.5),fliplr(1000*prctile(allV1nogo_rand_mean,2.5))],[1 0 0],'FaceAlpha',0.2,'EdgeAlpha',0);


end
title('V1 baseline: go (green) vs. nogo (red)')
xlabel('time (ms)');ylabel('average firing rate (Hz)');
xlim([-100 1000]);ylim([0 25]);
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
hold on;patch([ startpoint startpoint+dur3 startpoint+dur3 startpoint],[0 0 100 100],'b','FaceAlpha',0.1,'EdgeAlpha',0)
%hold on;patch([ startpoint startpoint+dur2 startpoint+dur2 startpoint],[0 0 100 100],'b','FaceAlpha',0.1,'EdgeAlpha',0)

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
figure;plot(timeoffset+((-Window):1*edgestep:(Window-1))/(samplingf/1000),1000*mean(allLMgo)/(edgestep/30),'Color',[0 0.8 0],'LineWidth',1.6);
hold on;plot(timeoffset+((-Window):1*edgestep:(Window-1))/(samplingf/1000),1000*mean(allLMnogo)/(edgestep/30),'Color',[1 0 0],'LineWidth',1.6);
title('LM baseline: go (green) vs. nogo (red)')
xlabel('time (ms)');ylabel('average firing rate (Hz)');
xlim([-100 1000]);ylim([0 25]);
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

%% imagesc : 2 figures: sorted off go or nogo responsivity
clm=0.5

[~,sorti_V1_go] = sort(sum(allV1go(:,floor(size(allV1go,2)/2):end),2) - sum(allV1go(:,1:floor(size(allV1go,2)/2)),2));
[~,sorti_V1_nogo] = sort(sum(allV1nogo(:,floor(size(allV1nogo,2)/2):end),2) - sum(allV1nogo(:,1:floor(size(allV1nogo,2)/2)),2));
[~,sorti_LM_go] = sort(sum(allLMgo(:,floor(size(allLMgo,2)/2):end),2) - sum(allLMgo(:,1:floor(size(allLMgo,2)/2)),2));
[~,sorti_LM_nogo] = sort(sum(allLMnogo(:,floor(size(allLMnogo,2)/2):end),2) - sum(allLMnogo(:,1:floor(size(allLMnogo,2)/2)),2));

startpoint=V1cells{1}.smb_centers.go(lag);
[~,sorti_LM_slc] = sort(sum(allLMgo(:,(round(size(allLMgo,2)/2)+startpoint/(edgestep/30)):(round(size(allLMgo,2)/2)+startpoint/(edgestep/30)+100/(edgestep/30))),2));
[~,sorti_LM_slc] = sort(sum(allLMnogo(:,(round(size(allLMnogo,2)/2)+startpoint/(edgestep/30)):(round(size(allLMnogo,2)/2)+startpoint/(edgestep/30)+100/(edgestep/30))),2));


%%%%%%% cells sorted based on go response size
indv_sorting = 2; % 2 is only for LM
sortingtype = 'go'; % only counts when indv_sorting = 0

f = figure;
f.Units = 'normalized';f.Position = [0.25 0.25 0.5 0.8];
if indv_sorting == 0
    f.Name = sprintf('V1 Cells sorted based on %s responses',sortingtype);
    s1=subplot(1,2,1);imagesc(timeoffset + ((-Window):1*edgestep:(Window-1))/(samplingf/1000),[],allV1go(eval(sprintf('sorti_V1_%s',sortingtype)),:));
elseif indv_sorting == 1 | indv_sorting == 2
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
s1.FontSize=8;
s1.XLabel.String = 'Time (ms)';
s1.YLabel.String = 'Cells';
hold on;patch([ startpoint startpoint+dur3 startpoint+dur3 startpoint],[s1.YLim(1) s1.YLim(1) s1.YLim(2) s1.YLim(2)],'b','FaceAlpha',0.1,'EdgeAlpha',0)

if indv_sorting == 0
    s2=subplot(1,2,2);imagesc(timeoffset + ((-Window):1*edgestep:(Window-1))/(samplingf/1000),[],allV1nogo(eval(sprintf('sorti_V1_%s',sortingtype)),:));
elseif indv_sorting == 1 | indv_sorting == 2
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
%hold on;patch([ startpoint startpoint+dur3 startpoint+dur3 startpoint],[s2.YLim(1) s2.YLim(1) s2.YLim(2) s2.YLim(2)],'b','FaceAlpha',0.1)
hold on;patch([ startpoint startpoint+dur3 startpoint+dur3 startpoint],[s2.YLim(1) s2.YLim(1) s2.YLim(2) s2.YLim(2)],'b','FaceAlpha',0.1,'EdgeAlpha',0)

f = figure;
f.Units = 'normalized';f.Position = [0.25 0.25 0.5 0.4];
if indv_sorting == 0
    f.Name = sprintf('V1 Cells sorted based on %s responses',sortingtype);
    s3=subplot(1,2,1);imagesc(timeoffset + ((-Window):1*edgestep:(Window-1))/(samplingf/1000),[],allLMgo(eval(sprintf('sorti_LM_%s',sortingtype)),:));
elseif indv_sorting == 1
    f.Name = 'V1 Cells sorted separately for different stimuli';
    s3=subplot(1,2,1);imagesc(timeoffset + ((-Window):1*edgestep:(Window-1))/(samplingf/1000),[],allLMgo(eval(sprintf('sorti_LM_%s','go')),:));
else
    f.Name = 'LM Cells sorted separately based on slc';
    s3=subplot(1,2,1);imagesc(timeoffset + ((-Window):1*edgestep:(Window-1))/(samplingf/1000),[],allLMgo(eval(sprintf('sorti_LM_%s','slc')),:));

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
hold on;patch([ startpoint startpoint+dur3 startpoint+dur3 startpoint],[s3.YLim(1) s3.YLim(1) s3.YLim(2) s3.YLim(2)],'b','FaceAlpha',0.1,'EdgeAlpha',0)

if indv_sorting == 0
    s4=subplot(1,2,2);imagesc(timeoffset + ((-Window):1*edgestep:(Window-1))/(samplingf/1000),[],allLMnogo(eval(sprintf('sorti_LM_%s',sortingtype)),:));
elseif indv_sorting == 1
    s4=subplot(1,2,2);imagesc(timeoffset + ((-Window):1*edgestep:(Window-1))/(samplingf/1000),[],allLMnogo(eval(sprintf('sorti_LM_%s','nogo')),:));
else
    s4=subplot(1,2,2);imagesc(timeoffset + ((-Window):1*edgestep:(Window-1))/(samplingf/1000),[],allLMnogo(eval(sprintf('sorti_LM_%s','slc')),:));

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
hold on;patch([ startpoint startpoint+dur3 startpoint+dur3 startpoint],[s4.YLim(1) s4.YLim(1) s4.YLim(2) s4.YLim(2)],'b','FaceAlpha',0.1,'EdgeAlpha',0)

set(gcf,'Color','w');
%% difference plots Version 1
dodifferenceplots=0;
if dodifferenceplots
V1_pool = [allV1nogo;allV1go];
LM_pool = [allLMnogo;allLMgo];

LMrand=nan(50000,size(allLMgo,2));
V1rand=nan(50000,size(allV1go,2));

for rrand = 1:50000
    selectV1=zeros(1,size(allV1go,1)+size(allV1nogo,1))';
    selectV1(randi(length(selectV1),[1,size(allV1nogo,1)]))=1;
    V1rand(rrand,:) = 1000*(mean(V1_pool((find(~selectV1)),:))-mean(V1_pool((find(selectV1)),:)))/(edgestep/30);
    
    selectLM=zeros(1,size(allLMgo,1)+size(allLMnogo,1))';
    selectLM(randi(length(selectLM),[1,size(allLMnogo,1)]))=1;
    LMrand(rrand,:) = 1000*(mean(LM_pool((find(~selectLM)),:))-mean(LM_pool((find(selectLM)),:)))/(edgestep/30);
end



figure;plot(((-Window):1*edgestep:(Window-1))/(samplingf/1000),1000*(mean(allV1go)-mean(allV1nogo))/(edgestep/30),'Color','b','LineWidth',1.6);
xpoints=[((-Window):1*edgestep:(Window-1))/(samplingf/1000) fliplr(((-Window):1*edgestep:(Window-1))/(samplingf/1000))];
ypoints=[prctile(V1rand,5),fliplr(prctile(V1rand,95))];
hold on;patch(xpoints,ypoints,'b','FaceAlpha',0.2,'EdgeAlpha',0)

title('V1 baseline: go - nogo')
xlabel('time (ms)');ylabel('averadodifferenceplotsge firing rate (Hz)');
xlim([-100 1000]);ylim([-10 10]);
axx = gca;
hold on; patch([0 500 500 0],[axx.YLim(1) axx.YLim(1) axx.YLim(2) axx.YLim(2)],[0.7 0.7 0.7],...
    'FaceAlpha',0.3,'EdgeColor','none')
axx.XTick=0:100:1000;
set(gcf,'Color','w')



figure;plot(((-Window):1*edgestep:(Window-1))/(samplingf/1000),1000*(mean(allLMgo)-mean(allLMnogo))/(edgestep/30),'Color','b','LineWidth',1.6);
xpoints=[((-Window):1*edgestep:(Window-1))/(samplingf/1000) fliplr(((-Window):1*edgestep:(Window-1))/(samplingf/1000))];
ypoints=[prctile(LMrand,5),fliplr(prctile(LMrand,95))];
hold on;patch(xpoints,ypoints,'b','FaceAlpha',0.2,'EdgeAlpha',0)

title('LM baseline: go - nogo')
xlabel('time (ms)');ylabel('average firing rate (Hz)');
xlim([-100 1000]);ylim([-10 10]);
axx = gca;
hold on; patch([0 500 500 0],[axx.YLim(1) axx.YLim(1) axx.YLim(2) axx.YLim(2)],[0.7 0.7 0.7],...
    'FaceAlpha',0.3,'EdgeColor','none')
axx.XTick=0:100:1000;
set(gcf,'Color','w')
%% difference plots - Version2
selindex=1;

LMrand=nan(1000,size(LMcells{1}.nbs.go,2));% 1000* time
V1rand=nan(1000,size(V1cells{1}.nbs.go,2));% 1000 * time
% for V1
for rrand = 1:1000  
    onerep_gonodif_percell_overtime_V1 = nan(numel(V1cells),size(V1cells{1}.nbs.go,2));
    for i= 1:numel(V1cells)
        
        V1_pool = [V1cells{i}.nbs.go;V1cells{i}.nbs.nogo];
        % LM_pool = [LMcells{i}.nbs.go;LMcells{i}.nbs.nogo];
        
        selectV1=zeros(1,size(V1cells{i}.nbs.go,1)+size(V1cells{i}.nbs.nogo,1))';
        selectV1(randi(length(selectV1),[1,size(V1cells{i}.nbs.nogo,1)]))=1;
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
    onerep_gonodif_percell_overtime_LM = nan(numel(LMcells),size(LMcells{1}.nbs.go,2));
    for i= 1:numel(LMcells)
        
        LM_pool = [LMcells{i}.nbs.go;LMcells{i}.nbs.nogo];
        % LM_pool = [LMcells{i}.nbs.go;LMcells{i}.nbs.nogo];
        
        selectLM=zeros(1,size(LMcells{i}.nbs.go,1)+size(LMcells{i}.nbs.nogo,1))';
        selectLM(randi(length(selectLM),[1,size(LMcells{i}.nbs.nogo,1)]))=1;
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