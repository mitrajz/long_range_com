% checking the delta time between different laser lags
% V1 and LM the same, go and nogo the same
ff=load(fullfile('/mnt/data/Mitra/figs/P2_L/bothdircombined',['cells_','FF','_pl20_an70_lw20_exG1_onlyC1_onlyS0_plstyle1.mat']));
fb=load(fullfile('/mnt/data/Mitra/figs/P2_L/bothdircombined',['cells_','FB','_pl20_an70_lw20_exG1_onlyC1_onlyS0_plstyle1.mat']));
figure;hist([cell2mat(cellfun(@(x) diff(x.smb_centers.go),ff.V1cells,'UniformOutput',0)),...
    cell2mat(cellfun(@(x) diff(x.smb_centers.go),fb.V1cells,'UniformOutput',0))])

%% IMPORTANT PARAMS
exptype = 'FB';
calcfbness = 1;
doshifted =0;
lsbs = 1; % only used when do shifted


%%%% load cells
load(fullfile('/mnt/data/Mitra/figs/P2_L/bothdircombined',['cells_',exptype,'_pl20_an65_lw20_exG1_onlyC1_onlyS0_plstyle1.mat']))
shifted = load(fullfile('/mnt/data/Mitra/figs/P2_L/bothdircombined',['cells_',exptype,'_pl20_an65-65_lw20_exG1_onlyC1_onlyS0_plstyle1.mat']));

cleanupcriteria.doclean = 1;
cleanupcriteria.significance = nan; % significance of activity instead of absolute value
cleanupcriteria.smb_activity = 2.5; % in Hz, 
cellremovecriteria.doclean = 1;
cellremovecriteria.stability = 1;
cellremovecriteria.lindriftTH = 0.5; % stability param - def:0.5
cellremovecriteria.activity = nan;
cellremovecriteria.responsiveness =[4 1];  % nan or several options:
cellremovecriteria.sigTH = 0.05; % responsiveness param - def:0.05
cellremovecriteria.minspkTH = 1; %responsiveness param (method 4 only)- def:1
cellremovecriteria.showexampleplots = 0; % order:targetgo, targetnogo, not_target go, not_targetnogo

cell_keep_ind = nan;
[LMcells,V1cells,Behcells,params] = cleanupcells(cleanupcriteria,LMcells,V1cells,Behcells,params);
[cell_keep_ind,LMcells,V1cells,Behcells,params] = removecells(cell_keep_ind,cellremovecriteria,LMcells,V1cells,Behcells,params);

% because cell_keep_ind is not set to nan, it uses remve_ind from above to
% remove shifted cells
[shifted.LMcells,shifted.V1cells,shifted.Behcells,shifted.params] = cleanupcells(cleanupcriteria,shifted.LMcells,shifted.V1cells,shifted.Behcells,shifted.params);
[cell_keep_ind,shifted.LMcells,shifted.V1cells,shifted.Behcells,shifted.params] = removecells(cell_keep_ind,cellremovecriteria,shifted.LMcells,shifted.V1cells,shifted.Behcells,shifted.params);



%%%%%% load lda file

ldafilename = [exptype,'lmodel_25Th.mat'];
filename = fullfile('/mnt/data/Mitra/figs/P2_L/bothdircombined/2020_02_24/LDA/version1',ldafilename);
load(filename);

allanimals = find(cellfun(@(x) isfield(x,'lmodel'),animalmodel))
numanimals =length(allanimals);
%%

if strcmp(exptype,'FF')
    targetcells = LMcells;
    shifted.targetcells = shifted.LMcells;
elseif strcmp(exptype,'FB')
    targetcells = V1cells;
    shifted.targetcells = shifted.V1cells;
end

timescalems = 65;
fo = fitoptions('Method','NonlinearLeastSquares','Lower',[-inf,-inf,0],'Upper',[inf,inf,inf],'StartPoint',[0.5,1,1000],'MaxIter',1000);
ftype = fittype('A*(exp(-x/tau)+B)','independent','x','options',fo);



avacr_go = [];
nfbness_go = [];
allcormat_go=[];
avacr_nogo = [];
nfbness_nogo = [];
allcormat_nogo=[];
%alltogether or indexed per animal?
fall = figure;
axall = gca;
rsall=[];
pall = [];

% example animal is 2: just change allanimals to 2
for animalnum=allanimals
    
    animal_avacr_go = [];
    animal_nfbness_go = [];
    
    animal_avacr_nogo = [];
    animal_nfbness_nogo = [];
    
    
    animalnum
    animalcellind=find(cellfun(@(x) x.simulcode,targetcells) == animalnum );
    if calcfbness % 
        fbness_go= (abs(nanmax(abs(animalmodel{animalnum}.lmodel.wt(:,:)))));% separate1:8 and 9:16?
        [~,fbnessgoind]=sort(fbness_go);
        fbness_nogo= (abs(nanmax(abs(animalmodel{animalnum}.lmodel.wt(:,:)))));% separate1:8 and 9:16?
        [~,fbnessnogoind]=sort(fbness_nogo);
    else
        fbness_go = nan(1,numel(animalcellind));
        fbness_nogo = nan(1,numel(animalcellind));
    end
    
    randeq = 0; % if zero, takes last. Seed is fixed  for random sampling
    targetcells(animalcellind) = equalizegonogotrials(targetcells(animalcellind),0,0,randeq);% frnormalize
    shifted.targetcells(animalcellind) = equalizegonogotrials(shifted.targetcells(animalcellind),0,0,randeq);% frnormalize
    
    for nnum=animalcellind

        if doshifted
            if lsbs == 1
                [cormat_go] = singleneuronACF_shifted(targetcells{nnum}.laAls.go,...
                    shifted.targetcells{nnum}.laAls.go);
                % just the first column of cormat_go has values: corresponding to lag 1
                avacr_go = [avacr_go; nanmean(cormat_go)];
            else
                [cormat_go] = singleneuronACF_shifted(targetcells{nnum}.laAbs.go,...
                    shifted.targetcells{nnum}.laAbs.go);
                % just the first column of cormat_go has values: corresponding to lag 1
                avacr_go = [avacr_go; nanmean(cormat_go)];
            end
        else
            [cormat_tc_go,cormat_go] = singleneuronACF(targetcells{nnum}.laAbs.go);
             allcormat_go = cat(3,allcormat_go,cormat_go);
            
            animal_avacr_go = [animal_avacr_go; nanmean(cormat_tc_go)];
            animal_nfbness_go = [animal_nfbness_go, fbness_go(find(animalcellind == nnum))];
            avacr_go = [avacr_go; nanmean(cormat_tc_go)];
            nfbness_go = [nfbness_go, fbness_go(find(animalcellind == nnum))];
        end
        
        
        
        %nogo
        if doshifted
            if lsbs == 1
                [cormat_nogo] = singleneuronACF_shifted(targetcells{nnum}.laAls.nogo,...
                    shifted.targetcells{nnum}.laAls.nogo);
                % just the first column of cormat_go has values: corresponding to lag 1
                avacr_nogo = [avacr_nogo; nanmean(cormat_nogo)];
            else
                [cormat_nogo] = singleneuronACF_shifted(targetcells{nnum}.laAbs.nogo,...
                    shifted.targetcells{nnum}.laAbs.nogo);
                % just the first column of cormat_go has values: corresponding to lag 1
                avacr_nogo = [avacr_nogo; nanmean(cormat_nogo)];
            end
        else
            [cormat_tc_nogo,cormat_nogo] = singleneuronACF(targetcells{nnum}.laAbs.nogo);
            allcormat_nogo = cat(3,allcormat_nogo,cormat_nogo);
            
            animal_avacr_nogo = [animal_avacr_nogo; nanmean(cormat_tc_nogo)];
            animal_nfbness_nogo = [animal_nfbness_nogo,fbness_nogo(find(animalcellind == nnum))];
            avacr_nogo = [avacr_nogo; nanmean(cormat_tc_nogo)];
            nfbness_nogo = [nfbness_nogo, fbness_nogo(find(animalcellind == nnum))];
        end
        %end
    end
    
    if ~doshifted
        [~,ind_go] = sort(animal_nfbness_go);
        [~,ind_nogo] = sort(animal_nfbness_nogo);
        mdl=fitlm([animal_nfbness_go(ind_go)],[animal_avacr_go((ind_go),2)']);
        x = [animal_nfbness_go(ind_go)];
        
        % uncomment if you want to see one example for one animal
       % scatter(x,[animal_avacr_go((ind_go),2)'],100,'k.')
        hold(axall,'on')
        xplot = 0:0.001:1;
        plot(xplot, ...
            mdl.Coefficients.Estimate(1)+ mdl.Coefficients.Estimate(2)*xplot,'k');
        [~,yp] = predict(mdl,xplot','Alpha',0.05,'Simultaneous',1);
        hold(axall,'on')
        patch([xplot,fliplr(xplot)],...
            [yp(:,1)',fliplr(yp(:,2)')],'k',...
            'EdgeAlpha',0,'FaceAlpha',0.1);
        hold(axall,'on')
        text(0.8,0.4+0.05*animalnum,['rsquared=',num2str(mdl.Rsquared.Ordinary)])
        text(0.5,0.4+0.05*animalnum,['p=',num2str(mdl.Coefficients.pValue(2))]) 
        %p-value for the t-statistic of the hypothesis test that the corresponding coefficient is equal to zero or not.
        hold(axall,'on')
        rsall(end+1) = (mdl.Rsquared.Ordinary);
        pall(end+1) = mdl.Coefficients.pValue(2); % this is p value of intercept
        

    end
end
ylim([-0.4,0.6])
text(0,0.6,['mean(ste) = ',num2str(mean(rsall)),'(',num2str(std(rsall)/sqrt(numel(rsall))),')'])

%% plots
% go vs nogo -- all cells
if doshifted
    figure;bar(1,nanmean(avacr_go(:,1)),0.1,'g','FaceAlpha',0.5,'EdgeAlpha',0);
    hold on; errorbar(1,nanmean(avacr_go(:,1)),2*nanstd(avacr_go(:,1))/sqrt(size(avacr_go(:,1),1)),'Color','g')
    hold on;bar(1.2,nanmean(avacr_nogo(:,1)),0.1,'r','FaceAlpha',0.5,'EdgeAlpha',0);
    hold on; errorbar(1.2,nanmean(avacr_nogo(:,1)),2*nanstd(avacr_nogo(:,1))/sqrt(size(avacr_nogo(:,1),1)),'Color','r')
    xlim([0.8,1.4])
    pv = ranksum(avacr_go(:,1),avacr_nogo(:,1));
    text(1,0.14,['ranksum p = ',num2str(pv)])
    title(['lsbs = ',num2str(lsbs)])
else
    xt = 0:7;
    figure;plot(xt(2:end),nanmean(avacr_go(:,2:end)),'g.');
    hold on; errorbar(xt(2:end),nanmean(avacr_go(:,2:end)),2*nanstd(avacr_go(:,2:end))/sqrt(size(avacr_go(:,2:end),1)),'Color','g')
    hold on;plot(xt(2:end),nanmean(avacr_nogo(:,2:end)),'r.');
    hold on; errorbar(xt(2:end),nanmean(avacr_nogo(:,2:end)),2*nanstd(avacr_nogo(:,2:end))/sqrt(size(avacr_nogo(:,2:end),1)),'Color','r')
    xlim([0 8]);
    ylim([-0.1 0.3])
    %%%%%%% model fits
    
    [goacmdl,gof] = autocorrelations_fit_helper(avacr_go,ftype,timescalems);
    figure;
    plot(timescalems*xt(2:end),goacmdl.A*(exp(-(1:7)/goacmdl.tau)+goacmdl.B),'-','Color','g')
    hold on;
    yp = predint(goacmdl,xt(2:end),0.95,'functional','on');
    patch(timescalems*[1:7,7:-1:1],...
        [yp(:,1)',fliplr(yp(:,2)')],'g',...
        'EdgeAlpha',0,'FaceAlpha',0.2);
    hold on;errorbar(timescalems*xt(2:end),nanmean(avacr_go(:,2:end)),2*nanstd(avacr_go(:,2:end))/sqrt(size(avacr_go(:,2:end),1)),...
        'LineStyle','none','Color','g')
    hold on;
    
    [nogoacmdl,gof] = autocorrelations_fit_helper(avacr_nogo,ftype,timescalems);
    plot(timescalems*xt(2:end),nogoacmdl.A*(exp(-(1:7)/nogoacmdl.tau)+nogoacmdl.B),'-','Color','r')
    hold on;
    yp = predint(nogoacmdl,xt(2:end),0.95,'functional','on');
    patch(timescalems*[1:7,7:-1:1],...
        [yp(:,1)',fliplr(yp(:,2)')],'r',...
        'EdgeAlpha',0,'FaceAlpha',0.2);
    hold on;errorbar(timescalems*xt(2:end),nanmean(avacr_nogo(:,2:end)),2*nanstd(avacr_nogo(:,2:end))/sqrt(size(avacr_nogo(:,2:end),1)),...
        'LineStyle','none','Color','r')
    
    
    %%%
    % go vs nogo -- separated based on fbness -- add go/nogo combined
    [~,ind_go] = sort(nfbness_go);
    [~,ind_nogo] = sort(nfbness_nogo);
    %%%
    mdl=fitlm([nfbness_go(ind_go),nfbness_nogo(ind_nogo)],[avacr_go((ind_go),2)',avacr_nogo((ind_nogo),2)']);
    x = [nfbness_go(ind_go),nfbness_nogo(ind_nogo)];
    figure;
    scatter(x,[avacr_go((ind_go),2)',avacr_nogo((ind_nogo),2)'],100,'k.')
    hold on;
    xplot = 0:0.001:1;
    plot(xplot, ...
        mdl.Coefficients.Estimate(1)+ mdl.Coefficients.Estimate(2)*xplot,'k');
    [~,yp] = predict(mdl,xplot','Alpha',0.05,'Simultaneous',1);
    patch([xplot,fliplr(xplot)],...
        [yp(:,1)',fliplr(yp(:,2)')],'k',...
        'EdgeAlpha',0,'FaceAlpha',0.2);
    text(0.8,0.4,['rsquared=',num2str(mdl.Rsquared.Ordinary)])
    
    %%%% go
    mdl=fitlm([nfbness_go(ind_go)],[avacr_go((ind_go),2)']);
    x = [nfbness_go(ind_go)];
    figure;
    scatter(x,[avacr_go((ind_go),2)'],100,'g.')
    hold on;
    xplot = 0:0.001:1;
    plot(xplot, ...
        mdl.Coefficients.Estimate(1)+ mdl.Coefficients.Estimate(2)*xplot,'g');
    [~,yp] = predict(mdl,xplot','Alpha',0.05,'Simultaneous',1);
    patch([xplot,fliplr(xplot)],...
        [yp(:,1)',fliplr(yp(:,2)')],'g',...
        'EdgeAlpha',0,'FaceAlpha',0.2);
    text(0.8,0.4,['rsquared=',num2str(mdl.Rsquared.Ordinary)])
    mdl=fitlm([nfbness_go(ind_go)],[avacr_go((ind_go),2)']);
    x = [nfbness_go(ind_go)];
    % same as above for nogo
    mdl=fitlm([nfbness_nogo(ind_nogo)],[avacr_nogo((ind_nogo),2)']);
    x = [nfbness_nogo(ind_nogo)];
    figure;
    scatter(x,[avacr_nogo((ind_nogo),2)'],100,'r.')
    hold on;
    xplot = 0:0.001:1;
    plot(xplot, ...
        mdl.Coefficients.Estimate(1)+ mdl.Coefficients.Estimate(2)*xplot,'r');
    [~,yp] = predict(mdl,xplot','Alpha',0.05,'Simultaneous',1);
    patch([xplot,fliplr(xplot)],...
        [yp(:,1)',fliplr(yp(:,2)')],'r',...
        'EdgeAlpha',0,'FaceAlpha',0.2);
    text(0.8,0.4,['rsquared=',num2str(mdl.Rsquared.Ordinary)])
end
%%
function [cormattc,cormat2] = singleneuronACF(ActivityVec)

min_n_trials = 10;

l=8;
cormat2 = nan(l,l);
if min(cellfun(@(x) numel(x),ActivityVec)) >= min_n_trials
    for i=1:l
        for j=1:l
            cv = cov(ActivityVec{i},ActivityVec{j})';
            v=sqrt(var(ActivityVec{i})*var(ActivityVec{j}));
            
            cormat2(i,j) = cv(1,2)/v(1,1);
        end
    end
end

% excluding corrlations between time points 1 and 2, that are ~50ms apart
% in some animals rather than 65

cormat2(1,2) = nan;
cormat2(2,1) = nan;

cormattc = nan(l,l);
currentrow = 0;
for ii=1:l
    currentrow = currentrow + 1;
    for jj =1:l
        if (jj-ii) >= 0
            cormattc(currentrow,jj-ii+1) = cormat2(ii,jj);
        end
    end
end
end
%% % this is only lag1 for main activity, and the shifted version of it
function [cormat2] = singleneuronACF_shifted(ActivityVec1,ActivityVec2)

min_n_trials = 10;
l=8;
cormat2 = nan(l,l);
if min(cellfun(@(x) numel(x),ActivityVec1)) >= min_n_trials
    for i=1:l
        for j=1
            cv = cov(ActivityVec1{i},ActivityVec2{i})';
            v=sqrt(var(ActivityVec1{i})*var(ActivityVec2{i}));
            
            cormat2(i,j) = cv(1,2)/v(1,1);
        end
    end
end

% excluding corrlations between time points 1 and 2, that are ~50ms apart
% in some animals rather than 65

cormat2(1,2) = nan;
cormat2(2,1) = nan;

end
%%
function [cormattc,cormat2] = pairwiseneuronCCF(ActivityVec1,ActivityVec2)

min_n_trials = 10;%14min(cellfun(@(x) numel(x),ActivityVec));

l=8;
cormat2 = nan(l,l);
if min(cellfun(@(x) numel(x),ActivityVec1)) >= min_n_trials
    for i=1:l
        for j=1:l
            cv = cov(ActivityVec1{i},ActivityVec2{j})';
            % v = corrcoef(ext(:,i),ext(:,i))';
            v=sqrt(var(ActivityVec1{i})*var(ActivityVec2{j}));
            
            %         cv = cov(ActivityVec{i}(1+end-min_n_trials:end),ActivityVec{j}(1+end-min_n_trials:end))';
            %         % v = corrcoef(ext(:,i),ext(:,i))';
            %        v=sqrt(var(ActivityVec{i}(1+end-min_n_trials:end))*var(ActivityVec{j}(1+end-min_n_trials:end)));
            cormat2(i,j) = cv(1,2)/v(1,1);
        end
    end
end

% excluding corrlations between time points 1 and 2, that are ~50ms apart
% in some animals rather than 65

cormat2(1,2) = nan;
cormat2(2,1) = nan;

cormattc = nan(l,l);
currentrow = 0;
for ii=1:l
    currentrow = currentrow + 1;
    for jj =1:l
        if (jj-ii) >= 0
            cormattc(currentrow,jj-ii+1) = cormat2(ii,jj);
        end
    end
end
end

%%
function [cormattc,cormat2] = singleneuronACF_laser(ActivityVec)
min_n_trials = 10;
l=8;
cormat2 = nan(l,l);
if min(cellfun(@(x) size(x,1),ActivityVec)) >= min_n_trials
    for i=1:l
        for j=i:min(i+3,l)
            cv = cov(ActivityVec{i}(:,1),ActivityVec{i}(:,j-i+1))';
            % v = corrcoef(ext(:,i),ext(:,i))';
            v=sqrt(var(ActivityVec{i}(:,1))*var(ActivityVec{i}(:,j-i+1)));
            
            %         cv = cov(ActivityVec{i}(1+end-min_n_trials:end),ActivityVec{j}(1+end-min_n_trials:end))';
            %         % v = corrcoef(ext(:,i),ext(:,i))';
            %        v=sqrt(var(ActivityVec{i}(1+end-min_n_trials:end))*var(ActivityVec{j}(1+end-min_n_trials:end)));
            cormat2(i,j) = cv(1,2)/v(1,1);
        end
    end
end

% excluding corrlations between time points 1 and 2, that are ~50ms apart
% in some animals rather than 65

cormat2(1,2) = nan;
cormat2(2,1) = nan;

cormattc = nan(l,l);
currentrow = 0;
for ii=1:l
    currentrow = currentrow + 1;
    for jj =1:l
        if (jj-ii) >= 0
            cormattc(currentrow,jj-ii+1) = cormat2(ii,jj);
        end
    end
end
end