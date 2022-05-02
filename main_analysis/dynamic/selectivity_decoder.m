% Decoding the stimulus identity from target activity in the presence and
% absence of long range input, using a lineardecoder(lda) with gamma.


Exptype = 'FB';

cleanupcriteria = struct;
cleanupcriteria.doclean = 1; % 1 or 0, perform any cleaning or not. If zero does nothing
cleanupcriteria.significance = nan; % significance of activity instead of absolute value
cleanupcriteria.smb_activity = 2.5; % in Hz

cellremovecriteria.doclean = 0; % 1 or 0, perform any cleaning or not. If zero does nothing
cell_keep_ind = nan; % important param: if given, indexes are no recalculated
cellremovecriteria.stability = 1;
cellremovecriteria.lindriftTH = 0.5; % stability param - def:0.5
cellremovecriteria.activity = nan;
cellremovecriteria.responsiveness =[4 1];  % nan or several options:
% first element: method
% second element: at least one of
% go or nogo, or both: 1 or 2
% def: 4,1
cellremovecriteria.sigTH = 0.05; % responsiveness param - def:0.05
cellremovecriteria.minspkTH = 1; %responsiveness param (method 4 only)- def:1
cellremovecriteria.showexampleplots = 0; % order:targetgo, targetnogo, not_target go, not_targetnogo

if strcmp(Exptype,'FB')
    load('cells_FB_pl20_an80_lw20_exG1_onlyC1_onlyS0_plstyle1.mat');
elseif strcmp(Exptype,'FF')
    load('cells_FF_pl20_an80_lw20_exG1_onlyC1_onlyS0_plstyle1.mat');
end

[LMcells,V1cells,Behcells,params] = cleanupcells(cleanupcriteria,LMcells,V1cells,Behcells,params);
[cell_keep_ind,LMcells,V1cells,Behcells,params] = removecells(cell_keep_ind,cellremovecriteria,LMcells,V1cells,Behcells,params);


%%% 
clc
timeWindowSec = 0.08;
if strcmp(Exptype,'FB')
    numanimals = 6;
elseif strcmp(Exptype,'FF')
    numanimals = 7;
end
nfold = 5; % 5 
nrep = 20; %20  
if strcmp(Exptype,'FB') && nfold == 5 && nrep == 20    
    best_gamma = [ 0.6000    0.5000    0.8000    0.9000    0.8000    0.7000];
elseif strcmp(Exptype,'FF')
     best_gamma = [ 0.9000    0.5000    1.0000    0.9000    0.6000    0.7000    0.8000];
end




dodif = 1;

numtimepoints = 8;
Mbs= nan(1,numtimepoints);
Mls= nan(1,numtimepoints);
Ebs= nan(1,numtimepoints);
Els= nan(1,numtimepoints);

bs_mat = nan(numanimals,numtimepoints);
ls_mat = nan(numanimals,numtimepoints);

A_bs_mat = nan(numanimals,numtimepoints);
A_ls_mat = nan(numanimals,numtimepoints);


for timepoint = 1:numtimepoints
    
    
    ClAc_bs = nan(numanimals,nrep);
    ClAc_ls = nan(numanimals,nrep);
    
    A_bs = nan(numanimals,nrep);
    A_ls = nan(numanimals,nrep);
    
    for animalnum=1:numanimals
        for rep=1:nrep
            
            if strcmp(Exptype,'FB')
                 targetcell = V1cells(find(cellfun(@(x) x.simulcode,V1cells) == animalnum ));
            elseif strcmp(Exptype,'FF')
                 targetcell = LMcells(find(cellfun(@(x) x.simulcode,LMcells) == animalnum ));
            end
           
            % ntrials*ncells:
            gobs = cell2mat(cellfun(@(x) x.laAbs.go{timepoint},targetcell,'UniformOutput',0));
            nogobs = cell2mat(cellfun(@(x) x.laAbs.nogo{timepoint},targetcell,'UniformOutput',0));
            gols = cell2mat(cellfun(@(x) x.laAls.go{timepoint},targetcell,'UniformOutput',0));
            nogols = cell2mat(cellfun(@(x) x.laAls.nogo{timepoint},targetcell,'UniformOutput',0));
            
            bsall = [gobs;nogobs];
            lsall = [gols;nogols];
            
            
            % equalize trials
            gobs(:,isnan(sum(bsall,1))) = [];
            nogobs(:,isnan(sum(bsall,1))) = [];
            gols(:,isnan(sum(lsall,1))) = [];
            nogols(:,isnan(sum(lsall,1))) = [];
            
            mintrialnum = min([size(gobs,1),size(nogobs,1),size(gols,1),size(nogols,1)]);
            gobs = gobs(1:mintrialnum,:);
            nogobs = nogobs(1:mintrialnum,:);
            
            %
            rng(rep)

            try

                cvind=crossvalind('Kfold',size(gobs,1)+size(nogobs,1),nfold);
                ClAc_bs_fold = [];
                ClAc_ls_fold = [];
                A_bs_fold = [];
                A_ls_fold = [];
                for fold=1:nfold
                    traininds = (find(cvind ~= fold));
                    testinds = (find(cvind == fold));
                    
                    PooledX = [gobs;nogobs];
                    PooledY = [zeros(size(gobs,1),1);ones(size(nogobs,1),1)];
                    mdl_bs = fitcdiscr(PooledX(traininds,:),...
                        PooledY(traininds,:),'DiscrimType','linear','Prior','uniform','Gamma',best_gamma(animalnum));
                    predicted_label = predict(mdl_bs,PooledX(testinds,:));
                    real_label = PooledY(testinds);
                    % check classification accuracy
                    A_bs_fold(end+1) = nanmean((nanmean(gobs,1)+nanmean(nogobs,1))/2);%sum(~xor(predicted_label,real_label))/numel(predicted_label);%nanmean((nanmean(gobs,1)+nanmean(nogobs,1))/2);%
                    ClAc_bs_fold(end+1) = sum(~xor(predicted_label,real_label))/numel(predicted_label);%nanmean((nanmean(gobs,1)+nanmean(nogobs,1))/2);%
                    
                    
                    predicted_label = predict(mdl_bs,[gols;nogols]);
                    real_label = [zeros(size(gols,1),1);ones(size(nogols,1),1)];
                    % check classification accuracy
                    A_ls_fold(end+1) = nanmean((nanmean(gols,1)+nanmean(nogols,1))/2);%sum(~xor(predicted_label,real_label))/numel(predicted_label);%nanmean((nanmean(gols,1)+nanmean(nogols,1))/2);%
                    ClAc_ls_fold(end+1) = sum(~xor(predicted_label,real_label))/numel(predicted_label);%nanmean((nanmean(gols,1)+nanmean(nogols,1))/2);%
                    
                end
                ClAc_bs(animalnum,rep) = nanmean(ClAc_bs_fold);
                ClAc_ls(animalnum,rep) = nanmean(ClAc_ls_fold);
                
                A_bs(animalnum,rep) = nanmean(A_bs_fold);
                A_ls(animalnum,rep) = nanmean(A_ls_fold);
            catch
            end
        end
        bs_mat(animalnum,timepoint) = nanmean(ClAc_bs(animalnum,:));
        ls_mat(animalnum,timepoint) = nanmean(ClAc_ls(animalnum,:));
        
        A_bs_mat(animalnum,timepoint) = nanmean(A_bs(animalnum,:));
        A_ls_mat(animalnum,timepoint) = nanmean(A_ls(animalnum,:));
    end
    
    [signrank(reshape(ClAc_bs,1,[]),reshape(ClAc_ls,1,[])),ranksum(reshape(ClAc_bs,1,[]),reshape(ClAc_ls,1,[]))]
    
    if dodif
        ClAc_bs = 100*(ClAc_ls - nanmean((ClAc_bs),2))./nanmean((ClAc_bs),2);
    end
    
    Mbs(timepoint) = nanmean(reshape(ClAc_bs,1,[]));
    Mls(timepoint) = nanmean(reshape(ClAc_ls,1,[]));
    
    
    Ebs(timepoint) = nanstd(reshape(ClAc_bs,1,[]))/sqrt(numel(ClAc_bs));
    Els(timepoint) = nanstd(reshape(ClAc_ls,1,[]))/sqrt(numel(ClAc_ls));
    
    if dodif
        Ebs(timepoint) = 2*Ebs(timepoint);
        Els(timepoint) = 2*Els(timepoint);
    end

end

if dodif
    figure;
    errorbar(1:8,Mbs,Ebs,'k.'); hold on
    line([0 8.5],[0 -0])
end


figure;shadedErrorBar([],Mbs,Ebs,'k',1);
 hold on; line([0 8.5],[0 -0])
%%
hold on;
shadedErrorBar([],Mls,Els,'b',1)

wdt = 0.2;
figure;
pb=bar((1:numtimepoints)-wdt/2,Mbs,wdt,'FaceColor','k','EdgeColor','none'); hold on
errorbar((1:numtimepoints)-wdt/2,Mbs,Ebs,'k.'); hold on
pl=bar((1:numtimepoints)+wdt/2,Mls,wdt,'FaceColor','b','EdgeColor','none'); hold on
errorbar((1:numtimepoints)+wdt/2,Mls,Els,'b.')

ylim([0.5,1])
ax1 = gca;
ax1.XTickLabel=[];
ylabel('stimulus decoding accurancy from V1 cells (Linear)');

xlim([0.5 3.5])
xlim([0 10])
ax1.XTick = [1 2 3];
ax1.XTickLabel={'0ms-70ms','60ms-130ms','120ms-190ms'};
legend([pb,pl],{'control','silencing'});
%%% silencing delta



%%% early/late
bs_early = reshape(bs_mat(:,1:2),1,[]);
bs_late = reshape(bs_mat(:,3:8),1,[]);
ls_early = reshape(ls_mat(:,1:2),1,[]);
ls_late = reshape(ls_mat(:,3:8),1,[]);

bs_early(find(isnan(bs_early))) = [];
bs_late(find(isnan(bs_late))) = [];
ls_early(find(isnan(ls_early))) = [];
ls_late(find(isnan(ls_late))) = [];

A_bs_early = reshape(A_bs_mat(:,1:2),1,[])/timeWindowSec;
A_bs_late = reshape(A_bs_mat(:,3:8),1,[])/timeWindowSec;
A_ls_early = reshape(A_ls_mat(:,1:2),1,[])/timeWindowSec;
A_ls_late = reshape(A_ls_mat(:,3:8),1,[])/timeWindowSec;

A_bs_early(find(isnan(A_bs_early))) = [];
A_bs_late(find(isnan(A_bs_late))) = [];
A_ls_early(find(isnan(A_ls_early))) = [];
A_ls_late(find(isnan(A_ls_late))) = [];

figure;
subplot(1,2,1)
hold on;bar(0.9,mean(bs_early),0.2,'FaceColor','k','FaceAlpha',0.5)
hold on;bar(1.1,mean(ls_early),0.2,'FaceColor','b','FaceAlpha',0.5)
hold on;bar(1.9,mean(bs_late),0.2,'FaceColor','k','FaceAlpha',0.5)
hold on;bar(2.1,mean(ls_late),0.2,'FaceColor','b','FaceAlpha',0.5)
hold on;line([.9,1.1], [mean(bs_early),mean(ls_early)]);
hold on;line([1.9,2.1], [mean(bs_late),mean(ls_late)]);
hold on;errorbar([.9,1.9],[mean(bs_early),mean(bs_late)],[std(bs_early)/sqrt(numel(bs_early)),std(bs_late)/sqrt(numel(bs_late))],'k')
hold on; errorbar([1.1,2.1],[mean(ls_early),mean(ls_late)],[std(ls_early)/sqrt(numel(ls_early)),std(ls_late)/sqrt(numel(ls_late))],'b')
text(1,.80,['signrank early bs ls = ',num2str(signrank(bs_early,ls_early))])
text(2,.85,['signrank late bs ls = ',num2str(signrank(bs_late,ls_late))])
ylim([0.5 1])
ylabel('decoding performance')


subplot(1,2,2)
hold on;bar(0.9,mean(A_bs_early),0.2,'FaceColor','k','FaceAlpha',0.5)
hold on;bar(1.1,mean(A_ls_early),0.2,'FaceColor','b','FaceAlpha',0.5)
hold on;bar(1.9,mean(A_bs_late),0.2,'FaceColor','k','FaceAlpha',0.5)
hold on;bar(2.1,mean(A_ls_late),0.2,'FaceColor','b','FaceAlpha',0.5)
hold on;line([.9,1.1], [mean(A_bs_early),mean(A_ls_early)]);
hold on;line([1.9,2.1], [mean(A_bs_late),mean(A_ls_late)]);
hold on;errorbar([.9,1.9],[mean(A_bs_early),mean(A_bs_late)],[std(A_bs_early)/sqrt(numel(A_bs_early)),std(A_bs_late)/sqrt(numel(A_bs_late))],'k')
hold on; errorbar([1.1,2.1],[mean(A_ls_early),mean(A_ls_late)],[std(A_ls_early)/sqrt(numel(A_ls_early)),std(A_ls_late)/sqrt(numel(A_ls_late))],'b')
text(1,.80,['signrank early bs ls = ',num2str(signrank(A_bs_early,A_ls_early))])
text(2,.85,['signrank late bs ls = ',num2str(signrank(A_bs_late,A_ls_late))])
ylabel('average fring rate, Hz')
