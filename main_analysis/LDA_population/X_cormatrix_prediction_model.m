% This script fits a glm to predict the FF_ij or FB_ij (ijth element of FB
% correlation atrix) from V1_ij or LM_ij, comibining all animals
% FFcormat
% FBcormat
% V1cormat
% LMcormat
% should make animalmodel files before running
%% Cleaning and removing from aggragate analysis:

%% load, clean up and cell remove params

normalizefr = 0; % 0 or 1
nfold = 10; % was 20
rngnum = 3; %
norm2obv = 1;

rootpath = '/mnt/data/Mitra/figs/P2_L/bothdircombined/2020_02_24/LDA/version2';
combinsizems = '150';
acbinsizems = '65';

nrep = 100;
Xstyle = 2;


% FF
exptype = 'FF';
[goV1linear, goLMlinear, goComcorlinear, nogoV1linear, nogoLMlinear, nogoComcorlinear]  = X_makelinearvars(exptype,...
    normalizefr,rootpath,combinsizems,acbinsizems,nrep,Xstyle);
[FF_Normalized_go_logl, FF_Normalized_nogo_logl,FF_go_beta,FF_nogo_beta] = X_makeCVLL(nfold,rngnum,norm2obv,...
    goV1linear, goLMlinear, goComcorlinear, nogoV1linear, nogoLMlinear, nogoComcorlinear);
% FB
exptype = 'FB';
[goV1linear, goLMlinear, goComcorlinear, nogoV1linear, nogoLMlinear, nogoComcorlinear]  = X_makelinearvars(exptype,...
    normalizefr,rootpath,combinsizems,acbinsizems,nrep,Xstyle);
[FB_Normalized_go_logl, FB_Normalized_nogo_logl,FB_go_beta,FB_nogo_beta] = X_makeCVLL(nfold,rngnum,norm2obv,...
    goV1linear, goLMlinear, goComcorlinear, nogoV1linear, nogoLMlinear, nogoComcorlinear);

%% plot mean and sem: TODO: % make units meaningful such as bits/etc

plotstyle = 2; % if 1: bar plot, if 2: strip plot
wdth = 0.2;
if plotstyle == 1
    
    % FF
    p1=bar(1-wdth,mean(FF_Normalized_go_logl),wdth,'FaceColor',[0 0.5 0],'EdgeColor','none');hold on;
    errorbar(1-wdth,mean(FF_Normalized_go_logl),std(FF_Normalized_go_logl)/sqrt(nfold),'.','Color',[0 0.5 0]);hold on
    
    p2=bar(2-wdth,mean(FF_Normalized_nogo_logl),wdth,'FaceColor',[0.5 0 0],'EdgeColor','none');hold on;
    errorbar(2-wdth,mean(FF_Normalized_nogo_logl),std(FF_Normalized_nogo_logl)/sqrt(nfold),'.','Color',[0.5 0 0]);hold on;
    % FB
    p3=bar(1+wdth,mean(FB_Normalized_go_logl),wdth,'FaceColor',[0 0.8 0],'EdgeColor','none');hold on;
    errorbar(1+wdth,mean(FB_Normalized_go_logl),std(FB_Normalized_go_logl)/sqrt(nfold),'.','Color',[0 0.8 0]);hold on
    
    p4=bar(2+wdth,mean(FB_Normalized_nogo_logl),wdth,'FaceColor',[0.8 0 0],'EdgeColor','none');hold on;
    errorbar(2+wdth,mean(FB_Normalized_nogo_logl),std(FB_Normalized_nogo_logl)/sqrt(nfold),'.','Color',[0.8 0 0]);hold on;
    
else
    markersize = 100;
    markerlpha = 0.3;
    figure
    scatter(1-wdth+(rand(1,length(FF_Normalized_go_logl))-0.5)*0.05,...
        FF_Normalized_go_logl,markersize,'.','MarkerEdgeColor',[0 0.5 0],'MarkerEdgeAlpha',markerlpha);hold on;
    p1=errorbar(1-wdth,mean(FF_Normalized_go_logl),2*std(FF_Normalized_go_logl)/sqrt(nfold),'.','Color',[0 0.5 0]);hold on
    bar(1-wdth,mean(FF_Normalized_go_logl),wdth,'FaceColor',[0 0.5 0],'FaceAlpha',0.3,'EdgeColor','none');hold on;
    
    scatter(1+wdth+(rand(1,length(FF_Normalized_nogo_logl))-0.5)*0.05,...
        FF_Normalized_nogo_logl,markersize,'.','MarkerEdgeColor',[0.5 0 0],'MarkerEdgeAlpha',markerlpha);hold on;
    p2=errorbar(1+wdth,mean(FF_Normalized_nogo_logl),2*std(FF_Normalized_nogo_logl)/sqrt(nfold),'.','Color',[0.5 0 0]);hold on;
    bar(1+wdth,mean(FF_Normalized_nogo_logl),wdth,'FaceColor',[0.5 0 0],'FaceAlpha',0.3,'EdgeColor','none');hold on;
    % FB
    scatter(2-wdth+(rand(1,length(FB_Normalized_go_logl))-0.5)*0.05,...
        FB_Normalized_go_logl,markersize,'.','MarkerEdgeColor',[0 0.8 0],'MarkerEdgeAlpha',markerlpha);hold on;
    p3=errorbar(2-wdth,mean(FB_Normalized_go_logl),2*std(FB_Normalized_go_logl)/sqrt(nfold),'.','Color',[0 0.8 0]);hold on
    bar(2-wdth,mean(FB_Normalized_go_logl),wdth,'FaceColor',[0 0.8 0],'FaceAlpha',0.3,'EdgeColor','none');hold on;
    
    scatter(2+wdth+(rand(1,length(FB_Normalized_nogo_logl))-0.5)*0.05,...
        FB_Normalized_nogo_logl,markersize,'.','MarkerEdgeColor',[0.8 0 0],'MarkerEdgeAlpha',markerlpha);hold on;
    p4=errorbar(2+wdth,mean(FB_Normalized_nogo_logl),2*std(FB_Normalized_nogo_logl)/sqrt(nfold),'.','Color',[0.8 0 0]);hold on;
    bar(2+wdth,mean(FB_Normalized_nogo_logl),wdth,'FaceColor',[0.8 0 0],'FaceAlpha',0.3,'EdgeColor','none');hold on;

end

legend([p1,p2,p3,p4],{'feedforward-go','feedforward-nogo','feedback-go','feedback-nogo'})

xlim([0.5 2.5])
ylabel('normalized log-likelihood / number of observations')
ax = gca;
ax.XTickLabel=[];
%% plotting coefficients:

fold = 1;
errorbarcoef = 2; % 1 = 1s.e. 2 = 2se (=95%ci)
%%%%%%%%%%%%%%%% source
figure
% FF
p1=errorbar(1-wdth,FF_go_beta.V1.estimate(fold),errorbarcoef*FF_go_beta.V1.se(fold),'.','Color',[0 0.5 0]);hold on
bar(1-wdth,FF_go_beta.V1.estimate(fold),wdth,'FaceColor',[0 0.5 0],'FaceAlpha',0.3,'EdgeColor','none');hold on;

p2=errorbar(1+wdth,FF_nogo_beta.V1.estimate(fold),errorbarcoef*FF_nogo_beta.V1.se(fold),'.','Color',[0.5 0 0]);hold on
bar(1+wdth,FF_nogo_beta.V1.estimate(fold),wdth,'FaceColor',[0.5 0 0],'FaceAlpha',0.3,'EdgeColor','none');hold on;
% FB
p3=errorbar(2-wdth,FB_go_beta.LM.estimate(fold),errorbarcoef*FB_go_beta.LM.se(fold),'.','Color',[0 0.8 0]);hold on
bar(2-wdth,FB_go_beta.LM.estimate(fold),wdth,'FaceColor',[0 0.8 0],'FaceAlpha',0.3,'EdgeColor','none');hold on;

p4=errorbar(2+wdth,FB_nogo_beta.LM.estimate(fold),errorbarcoef*FB_nogo_beta.LM.se(fold),'.','Color',[0.8 0 0]);hold on
bar(2+wdth,FB_nogo_beta.LM.estimate(fold),wdth,'FaceColor',[0.8 0 0],'FaceAlpha',0.3,'EdgeColor','none');hold on;
title('source')
ylim([-1 2.5])
text([1-wdth 1+wdth 2-wdth 2+wdth],[2 1.8 2 1.8],...
  arrayfun(@(x) num2str(x),...
  [ FF_go_beta.V1.p(fold),FF_nogo_beta.V1.p(fold),FB_go_beta.LM.p(fold),FB_nogo_beta.LM.p(fold)],'UniformOutput',0) )


%%%%%%%%%%%%%%%% target
figure
% FF
p1=errorbar(1-wdth,FF_go_beta.LM.estimate(fold),errorbarcoef*FF_go_beta.LM.se(fold),'.','Color',[0 0.5 0]);hold on
bar(1-wdth,FF_go_beta.LM.estimate(fold),wdth,'FaceColor',[0 0.5 0],'FaceAlpha',0.3,'EdgeColor','none');hold on;

p2=errorbar(1+wdth,FF_nogo_beta.LM.estimate(fold),errorbarcoef*FF_nogo_beta.LM.se(fold),'.','Color',[0.5 0 0]);hold on
bar(1+wdth,FF_nogo_beta.LM.estimate(fold),wdth,'FaceColor',[0.5 0 0],'FaceAlpha',0.3,'EdgeColor','none');hold on;
% FB
p3=errorbar(2-wdth,FB_go_beta.V1.estimate(fold),errorbarcoef*FB_go_beta.V1.se(fold),'.','Color',[0 0.8 0]);hold on
bar(2-wdth,FB_go_beta.V1.estimate(fold),wdth,'FaceColor',[0 0.8 0],'FaceAlpha',0.3,'EdgeColor','none');hold on;

p4=errorbar(2+wdth,FB_nogo_beta.V1.estimate(fold),errorbarcoef*FB_nogo_beta.V1.se(fold),'.','Color',[0.8 0 0]);hold on
bar(2+wdth,FB_nogo_beta.V1.estimate(fold),wdth,'FaceColor',[0.8 0 0],'FaceAlpha',0.3,'EdgeColor','none');hold on;
title('target')
ylim([-1 2.5])
text([1-wdth 1+wdth 2-wdth 2+wdth],[2 1.8 2 1.8],...
  arrayfun(@(x) num2str(x),...
  [ FF_go_beta.LM.p(fold),FF_nogo_beta.LM.p(fold),FB_go_beta.V1.p(fold),FB_nogo_beta.V1.p(fold)],'UniformOutput',0) )


% bias
figure
% FF
p1=errorbar(1-wdth,FF_go_beta.bias.estimate(fold),errorbarcoef*FF_go_beta.bias.se(fold),'.','Color',[0 0.5 0]);hold on
bar(1-wdth,FF_go_beta.bias.estimate(fold),wdth,'FaceColor',[0 0.5 0],'FaceAlpha',0.3,'EdgeColor','none');hold on;

p2=errorbar(1+wdth,FF_nogo_beta.bias.estimate(fold),errorbarcoef*FF_nogo_beta.bias.se(fold),'.','Color',[0.5 0 0]);hold on
bar(1+wdth,FF_nogo_beta.bias.estimate(fold),wdth,'FaceColor',[0.5 0 0],'FaceAlpha',0.3,'EdgeColor','none');hold on;
% FB
p3=errorbar(2-wdth,FB_go_beta.bias.estimate(fold),errorbarcoef*FB_go_beta.bias.se(fold),'.','Color',[0 0.8 0]);hold on
bar(2-wdth,FB_go_beta.bias.estimate(fold),wdth,'FaceColor',[0 0.8 0],'FaceAlpha',0.3,'EdgeColor','none');hold on;

p4=errorbar(2+wdth,FB_nogo_beta.bias.estimate(fold),errorbarcoef*FB_nogo_beta.bias.se(fold),'.','Color',[0.8 0 0]);hold on
bar(2+wdth,FB_nogo_beta.bias.estimate(fold),wdth,'FaceColor',[0.8 0 0],'FaceAlpha',0.3,'EdgeColor','none');hold on;
title('bias')

ylim([-2 2])
text([1-wdth 1+wdth 2-wdth 2+wdth],[1.5 1.3 1.5 1.3],...
  arrayfun(@(x) num2str(x),...
  [ FF_go_beta.bias.p(fold),FF_nogo_beta.bias.p(fold),FB_go_beta.bias.p(fold),FB_nogo_beta.bias.p(fold)],'UniformOutput',0) )


%%
function [Normalized_go_logl, Normalized_nogo_logl,go_beta,nogo_beta] = X_makeCVLL(nfold,rngnum,norm2obv,...
    goV1linear, goLMlinear, goComcorlinear, nogoV1linear, nogoLMlinear, nogoComcorlinear)
%%% crossvalidated loglikelihood (stratified) normalized to null model, per observation
% model is glm with a linear link. Log likelihood calculated for a gaussian
% distribution


rng(rngnum)
cvind=crossvalind('Kfold',size(goV1linear,2),nfold);

go_beta.V1.estimate = nan(1,nfold);
nogo_beta.V1.estimate = nan(1,nfold);
go_beta.LM.estimate = nan(1,nfold);
nogo_beta.LM.estimate = nan(1,nfold);
go_beta.bias.estimate = nan(1,nfold);
nogo_beta.bias.estimate = nan(1,nfold);

go_beta.V1.se = nan(1,nfold);
nogo_beta.V1.se = nan(1,nfold);
go_beta.LM.se = nan(1,nfold);
nogo_beta.LM.se = nan(1,nfold);
go_beta.bias.se = nan(1,nfold);
nogo_beta.bias.se = nan(1,nfold);

go_beta.V1.p = nan(1,nfold);
nogo_beta.V1.p = nan(1,nfold);
go_beta.LM.p = nan(1,nfold);
nogo_beta.LM.p = nan(1,nfold);
go_beta.bias.p = nan(1,nfold);
nogo_beta.bias.p = nan(1,nfold);


go_logl = nan(1,nfold);
go_logl_null = nan(1,nfold);
Normalized_go_logl = nan(1,nfold);
nogo_logl = nan(1,nfold);
nogo_logl_null = nan(1,nfold);
Normalized_nogo_logl = nan(1,nfold);

for fold=1:nfold
    traininds = (find(cvind ~= fold));
    testinds = (find(cvind == fold));
    
    %%%%%%%%%%%%% go
    % make main model and null model from training set
    gomdl = fitglm([reshape(goV1linear(:,traininds),1,[])',reshape(goLMlinear(:,traininds),1,[])'],...
        reshape(goComcorlinear(:,traininds),1,[])');
    
    go_beta.bias.estimate(fold) = gomdl.Coefficients.Estimate(1);
    go_beta.V1.estimate(fold)  = gomdl.Coefficients.Estimate(2);
    go_beta.LM.estimate(fold)  = gomdl.Coefficients.Estimate(3);
    go_beta.bias.se(fold)  = gomdl.Coefficients.SE(1);
    go_beta.V1.se(fold)  = gomdl.Coefficients.SE(2);
    go_beta.LM.se(fold)  = gomdl.Coefficients.SE(3);
    go_beta.bias.p(fold)  = gomdl.Coefficients.pValue(1);
    go_beta.V1.p(fold)  = gomdl.Coefficients.pValue(2);
    go_beta.LM.p(fold)  = gomdl.Coefficients.pValue(3);
    
    
    
   
    testx = [reshape(goV1linear(:,testinds),1,[])',reshape(goLMlinear(:,testinds),1,[])'];
    testy = reshape(goComcorlinear(:,testinds),1,[])';

    sigma = sqrt(gomdl.Dispersion);
    yhat = predict(gomdl,testx);
    if numel(find(isnan(yhat)))
        error
    end
    yhatnull = repmat(gomdl.Coefficients.Estimate(1),numel(yhat),1);
    sigmanull = std(reshape(goComcorlinear(:,traininds),1,[])');
    N = length(testy);

    go_logl(fold) = (-N/2)*log(2*pi*sigma^2) - ((1/(2*sigma^2)) * sum((yhat - testy).^2));
    go_logl_null(fold) = (-N/2)*log(2*pi*sigmanull^2) - ((1/(2*sigmanull^2)) * sum((yhatnull - testy).^2));

    if norm2obv
        Normalized_go_logl(fold) = (go_logl(fold) - go_logl_null(fold))/N;
    else
        Normalized_go_logl(fold) = (go_logl(fold) - go_logl_null(fold));
    end
    
    %%%%%%%%%%%%% nogo
    % make main model and null model from training set
    nogomdl = fitglm([reshape(nogoV1linear(:,traininds),1,[])',reshape(nogoLMlinear(:,traininds),1,[])'],...
        reshape(nogoComcorlinear(:,traininds),1,[])');
    
    nogo_beta.bias.estimate(fold)  = nogomdl.Coefficients.Estimate(1);
    nogo_beta.V1.estimate(fold)  = nogomdl.Coefficients.Estimate(2);
    nogo_beta.LM.estimate(fold)  = nogomdl.Coefficients.Estimate(3);
    nogo_beta.bias.se(fold)  = nogomdl.Coefficients.SE(1);
    nogo_beta.V1.se(fold)  = nogomdl.Coefficients.SE(2);
    nogo_beta.LM.se(fold)  = nogomdl.Coefficients.SE(3);
    nogo_beta.bias.p(fold)  = nogomdl.Coefficients.pValue(1);
    nogo_beta.V1.p(fold)  = nogomdl.Coefficients.pValue(2);
    nogo_beta.LM.p(fold)  = nogomdl.Coefficients.pValue(3);
    
  
    % log likelihood of test set under main model and test model
    testx = [reshape(nogoV1linear(:,testinds),1,[])',reshape(nogoLMlinear(:,testinds),1,[])'];
    testy = reshape(nogoComcorlinear(:,testinds),1,[])';
    sigma = sqrt(nogomdl.Dispersion);
    yhat = predict(nogomdl,testx);
    if numel(find(isnan(yhat)))
        error
    end
    yhatnull = repmat(nogomdl.Coefficients.Estimate(1),numel(yhat),1);
    % calculating sigma null like bellow is equivalent to training a model
    % with zero predictor and calcultaing its sigma:
    sigmanull = std(reshape(nogoComcorlinear(:,traininds),1,[])');
    N = length(testy);
    

    nogo_logl(fold) = (-N/2)*log(2*pi*sigma^2) - ((1/(2*sigma^2)) * sum((yhat - testy).^2));
    nogo_logl_null(fold) = (-N/2)*log(2*pi*sigmanull^2) - ((1/(2*sigmanull^2)) * sum((yhatnull - testy).^2));

    if norm2obv
        Normalized_nogo_logl(fold) = (nogo_logl(fold) - nogo_logl_null(fold))/N;
    else
        Normalized_nogo_logl(fold) = (nogo_logl(fold) - nogo_logl_null(fold));
    end
end


%%% normalize to get meaningful units
% nobsv =  size(goV1linear,2)*size(goV1linear,1)/nfold; % approximate number of observations per fold;
% if norm2obv
%     Normalized_go_logl = Normalized_go_logl/nobsv;
%     Normalized_nogo_logl = Normalized_nogo_logl/nobsv;
% end
end

%%
function [goV1linear, goLMlinear, goComcorlinear, nogoV1linear, nogoLMlinear, nogoComcorlinear]  = ...
    X_makelinearvars(exptype,normalization,rootpath,combinsizems,acbinsizems,nrep,Xstyle)

% no need to load cells. animalcells contain all
% used for activity only
ldafilename = [exptype,'lmodel_25Th_',acbinsizems,'ms_X_style',num2str(Xstyle),'_nrep',num2str(nrep),'.mat'];
fulfilename = fullfile(rootpath,ldafilename);
load(fulfilename)

allanimals = find(cellfun(@(x) isfield(x,'lmodel'),animalmodel))
numanimals =length(allanimals);


% used for com only
ldafilename = [exptype,'lmodel_25Th_',combinsizems,'ms_X_style',num2str(Xstyle),'_nrep',num2str(nrep),'.mat'];
fulfilename = fullfile(rootpath,ldafilename);
B=load(fulfilename);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% targetcell and notargetcell should have been specified before here
% load, clean and remove before here (aggregate analysis)

clear  goV1linear goLMlinear goComcorlinear nogoV1linear nogoLMlinear nogoComcorlinear
for animalnum = 1:numanimals
 
    % Here bs go and bs nogo have the same number of trials. (same as ls go
    % and ls nogo)
    if normalization == 1
        if strcmp(exptype,'FF')
            accormat_V1 =animalmodel{animalnum}.lmodel.acmatgo_notarget_n1;
            accormat_LM = animalmodel{animalnum}.lmodel.acmatgo_n1;
        elseif strcmp(exptype,'FB')
            accormat_V1 =animalmodel{animalnum}.lmodel.acmatgo_n1;
            accormat_LM = animalmodel{animalnum}.lmodel.acmatgo_notarget_n1;
        end
    elseif normalization == 0
        if strcmp(exptype,'FF')
            accormat_V1 =animalmodel{animalnum}.lmodel.acmatgo_notarget_n0;
            accormat_LM = animalmodel{animalnum}.lmodel.acmatgo_n0;
        elseif strcmp(exptype,'FB')
            accormat_V1 =animalmodel{animalnum}.lmodel.acmatgo_n0;
            accormat_LM = animalmodel{animalnum}.lmodel.acmatgo_notarget_n0;
        end
    end
    Comcormat = B.animalmodel{animalnum}.lmodel.Cmatgo_lda;
    % linearize the 3 maps: 1 row for each animal
    goV1linear(animalnum,:) = [accormat_V1(1,2:8),accormat_V1(2,3:8),accormat_V1(3,4:8),accormat_V1(4,5:8),...
        accormat_V1(5,6:8),accormat_V1(6,7:8),accormat_V1(7,8)];
    goLMlinear(animalnum,:) = [accormat_LM(1,2:8),accormat_LM(2,3:8),accormat_LM(3,4:8),accormat_LM(4,5:8),...
        accormat_LM(5,6:8),accormat_LM(6,7:8),accormat_LM(7,8)];
    goComcorlinear(animalnum,:) = [Comcormat(1,2:8),Comcormat(2,3:8),Comcormat(3,4:8),Comcormat(4,5:8),Comcormat(5,6:8),Comcormat(6,7:8),Comcormat(7,8)];
    
    
    % nogo: 
    if normalization == 1
        if strcmp(exptype,'FF')
            accormat_V1 =animalmodel{animalnum}.lmodel.acmatnogo_notarget_n1;
            accormat_LM = animalmodel{animalnum}.lmodel.acmatnogo_n1;
        elseif strcmp(exptype,'FB')
            accormat_V1 =animalmodel{animalnum}.lmodel.acmatnogo_n1;
            accormat_LM = animalmodel{animalnum}.lmodel.acmatnogo_notarget_n1;
        end
    elseif normalization == 0
        if strcmp(exptype,'FF')
            accormat_V1 =animalmodel{animalnum}.lmodel.acmatnogo_notarget_n0;
            accormat_LM = animalmodel{animalnum}.lmodel.acmatnogo_n0;
        elseif strcmp(exptype,'FB')
            accormat_V1 =animalmodel{animalnum}.lmodel.acmatnogo_n0;
            accormat_LM = animalmodel{animalnum}.lmodel.acmatnogo_notarget_n0;
        end
    end
    Comcormat = B.animalmodel{animalnum}.lmodel.Cmatnogo_lda;

    % linearize the 3 maps: 1 row for each animal
    nogoV1linear(animalnum,:) = [accormat_V1(1,2:8),accormat_V1(2,3:8),accormat_V1(3,4:8),accormat_V1(4,5:8),...
        accormat_V1(5,6:8),accormat_V1(6,7:8),accormat_V1(7,8)];
    nogoLMlinear(animalnum,:) = [accormat_LM(1,2:8),accormat_LM(2,3:8),accormat_LM(3,4:8),accormat_LM(4,5:8),...
        accormat_LM(5,6:8),accormat_LM(6,7:8),accormat_LM(7,8)];    
    nogoComcorlinear(animalnum,:) = [Comcormat(1,2:8),Comcormat(2,3:8),Comcormat(3,4:8),Comcormat(4,5:8),Comcormat(5,6:8),Comcormat(6,7:8),Comcormat(7,8)];
    %
    
end
nananimals = find(isnan(sum(nogoComcorlinear,2)+sum(goComcorlinear,2)));
goV1linear(nananimals,:)=[];
goLMlinear(nananimals,:)=[];
goComcorlinear(nananimals,:)=[];
nogoV1linear(nananimals,:)=[];
nogoLMlinear(nananimals,:)=[];
nogoComcorlinear(nananimals,:)=[];
end
