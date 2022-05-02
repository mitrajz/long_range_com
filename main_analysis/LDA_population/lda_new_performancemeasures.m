% distribution of error per animals
% This is just checking that lda is not purely noise. It does not show the
% shuffle test/ consistency of angles over bootstrapping. That needs to be
% done separately
% questions: - peranimal does lda perform better than chance?
%            - systematic difference in dv error between go and nogo?
%            - systematic difference in cv between different lags?

%% preprocessing: load lda file and see which animals have models
rootdir = '/mnt/data/Mitra/figs/P2_L/bothdircombined/2020_02_24/LDA/version1';
filename = 'FBlmodel_25Th.mat';
load(fullfile(rootdir,filename));
allanimals = find(cellfun(@(x) isfield(x,'lmodel'),animalmodel))

addshuffle = 1;
exptype = 'FF';
if addshuffle
    shf = load(fullfile('/mnt/data/Mitra/figs/P2_L/bothdircombined/2020_02_24/LDA/version1',[exptype,'lmodel_25Th_shuffle.mat']));
end

%% 1 - averaged over all lags for each animal: look at the spread over folds
style = 2;% DEFAULT:2 1: boxplot 2: mean+2*ste 95%CI
separateg_ng = 0; % DEFAULT: do both if 1: the number on top of each animal is signrank(go,nogo) - 0 or 1
figure;
for animal = allanimals
    if separateg_ng
        animalvec1 = nanmean(cell2mat(animalmodel{animal}.lmodel.Tefolds_go'),1);
        animalvec2 = nanmean(cell2mat(animalmodel{animal}.lmodel.Tefolds_nogo'),1)';
        if style == 1
            hold on; boxplot(animalvec1,...
                'PlotStyle','compact','Positions',animal+zeros(1,length(animalvec1)),'Colors','g','MedianStyle','line');
            hold on; boxplot(animalvec2,...
                'PlotStyle','compact','Positions',animal+0.2+zeros(1,length(animalvec2)),'Colors','r','MedianStyle','line');
            
        elseif style == 2
            hold on; errorbar(animal,nanmean(animalvec1),2*nanstd(animalvec1)/sqrt(length(animalvec1)),...
                'Color','g')
            hold on;scatter(animal,nanmean(animalvec1),100,'g.')
            hold on; errorbar(animal+0.2,nanmean(animalvec2),2*nanstd(animalvec2)/sqrt(length(animalvec2)),...
                'Color','r')
            hold on;scatter(animal+0.2,nanmean(animalvec2),100,'r.')
            text(animal,0.8,num2str(signrank(animalvec2,animalvec1)))
        end
    else
        animalvec = nanmean([nanmean(cell2mat(animalmodel{animal}.lmodel.Tefolds_go'),1)',...
            nanmean(cell2mat(animalmodel{animal}.lmodel.Tefolds_nogo'),1)'],2);
        if style == 1
            % average over folds: with median hyperparams (or as is saved in lmodel)
            hold on; boxplot(animalvec,...
                'PlotStyle','compact','Positions',animal+zeros(1,length(animalvec)),'Colors','g','MedianStyle','line');
        elseif style == 2
            hold on; errorbar(animal,nanmean(animalvec),2*nanstd(animalvec)/sqrt(length(animalvec)),...
                'Color','k')
            hold on;scatter(animal,nanmean(animalvec),100,'k.')
            
        end
    end
    
    
    if addshuffle
        % For each animal, for a given condition (go/nogo), there are 8
        % lags. For each lag there are 100 shuffle models, and for each
        % shuffle model, 10 folds. The error across the 10 folds are
        % averaged, and the distribution of average per fold error across
        % 100 shuffles are plotted, averaged for go/nogo. The error bars
        % are mean+-2*std  

        % size: (nboot)*nlags
        goers = cell2mat(cellfun(@(y) cellfun(@(x) nanmean(x),y.cv.test_er), shf.animalmodel{animal}.lmodel.shuffle_go,'UniformOutput',0)');
        nogoers = cell2mat(cellfun(@(y) cellfun(@(x) nanmean(x),y.cv.test_er), shf.animalmodel{animal}.lmodel.shuffle_nogo,'UniformOutput',0)');
        % size: (nboot)*1 (average over lags)
        sh_animalvec=nanmean([nanmean(goers,1);nanmean(nogoers,1)],1);
        
        hold on; errorbar(animal,nanmean(sh_animalvec),2*nanstd(sh_animalvec)/1,...
            'Color',[0.5 0.5 0.5])
        hold on;scatter(animal,nanmean(sh_animalvec),100,'.','MarkerEdgeColor',[0.5 0.5 0.5])
        
    end
end
xlim([0,8])
hold on;line([0 10],[0.5 0.5],'Color','k','LineStyle','--')
ylim([0 1])
%% 2- systemetic differences between go and nogo (can also be seen from above)
% per animal differences are characterized above. Here is over all animals,
% average test error
figure;
allgomeans = [];
allnogomeans = [];
for animal = allanimals
    animalvec1 = nanmean(cell2mat(animalmodel{animal}.lmodel.Tefolds_go'),1);
    animalvec2 = nanmean(cell2mat(animalmodel{animal}.lmodel.Tefolds_nogo'),1)';
    allgomeans(end+1) = nanmean(animalvec1);
    allnogomeans(end+1) = nanmean(animalvec2);
    hold on; line([0 1],[nanmean(animalvec1),nanmean(animalvec2)],'Color','k','LineStyle','--')
    hold on;
    scatter(0,nanmean(animalvec1),100,'k.');
    hold on;
    scatter(1,nanmean(animalvec2),100,'k.')
end
xlim([-1 2])
ylim([0 1])
hold on; text(0,0.7,['signrank p = ',num2str(signrank(allgomeans,allnogomeans))])
%% 3- systemetic differences between lags
%  over all animals, average test error

separateg_ng = 0; % DO BOTH

figure;
allgomeans = [];
allnogomeans = [];
allanimalvec = nan(max(allanimals),8);

for lag=1:8
    % averaging over animals and folds
    animalvec1 = cellfun(@(x) nanmean(x.lmodel.Tefolds_go{lag}),animalmodel(allanimals));
    animalvec2 = cellfun(@(x) nanmean(x.lmodel.Tefolds_nogo{lag}),animalmodel(allanimals));
    % averaging over go and nogo
    animalvec = nanmean([animalvec1;animalvec2],1);
    allanimalvec(:,lag) = animalvec;
    if separateg_ng
        hold on; errorbar(lag,nanmean(animalvec1),2*nanstd(animalvec1)/sqrt(length(animalvec1)),...
            'Color','g')
        hold on;scatter(lag,nanmean(animalvec1),100,'g.')
        hold on; errorbar(lag+0.2,nanmean(animalvec2),2*nanstd(animalvec2)/sqrt(length(animalvec2)),...
            'Color','r')
        hold on;scatter(lag+0.2,nanmean(animalvec2),100,'r.')
        text(lag,0.8,num2str(signrank(animalvec2,animalvec1)))
    else
        
        hold on; errorbar(lag,nanmean(animalvec),2*nanstd(animalvec)/sqrt(length(animalvec)),...
            'Color','k')
        hold on;scatter(lag,nanmean(animalvec),100,'k.')
        
    end
end

    xlim([-1 10])
    ylim([0 1])
    hold on
    text(0,0.7,['N in each group = ',num2str(length(allanimals))])
p = anova1(allanimalvec)