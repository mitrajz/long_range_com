
exptype = 'FB'; 

load(['/mnt/data/Mitra/figs/P2_L/bothdircombined/2020_02_24/LDA/version2/',exptype,'lmodel_25Th_150ms_X_style0_nrep100.mat'])
sh=load(['/mnt/data/Mitra/figs/P2_L/bothdircombined/2020_02_24/LDA/version2/shuffle_',exptype,'lmodel_25Th_150ms_X_style0_nrep100.mat']);
%%
allanimals = find(cellfun(@(x) isfield(x,'lmodel'),animalmodel))
numanimals =length(allanimals);

govec =nan(numanimals,8);
nogovec=nan(numanimals,8);
shgovec =nan(numanimals,8,100,2);
shnogovec=nan(numanimals,8,100,2);


add_chance_ci = 0;
distMeasure = 'batta'; % batta, lda-norm, mu-norm
for animalnum= allanimals
    
    for i=1:8
        if strcmp(distMeasure,'batta')
            govec(animalnum,i) = nanmean((cellfun(@(x) x.part{1}.go.batta_dist{i} + x.part{2}.go.batta_dist{i},animalmodel{animalnum}.lmodel.rep))/2);
            nogovec(animalnum,i) = nanmean((cellfun(@(x) x.part{1}.nogo.batta_dist{i} + x.part{2}.nogo.batta_dist{i},animalmodel{animalnum}.lmodel.rep))/2);
            shgovec(animalnum,i,:,:) = [(cellfun(@(x) x.part{1}.go.batta_dist{i} ,sh.animalmodel{animalnum}.lmodel.rep));(cellfun(@(x) x.part{2}.go.batta_dist{i} ,sh.animalmodel{animalnum}.lmodel.rep))]';
            shnogovec(animalnum,i,:,:) = [(cellfun(@(x) x.part{1}.nogo.batta_dist{i} ,sh.animalmodel{animalnum}.lmodel.rep));(cellfun(@(x) x.part{2}.nogo.batta_dist{i} ,sh.animalmodel{animalnum}.lmodel.rep))]';
            
        elseif strcmp(distMeasure,'lda-norm')
            govec(animalnum,i) = nanmean((cellfun(@(x) x.part{1}.go.lda_norm_dist{i} + x.part{2}.go.lda_norm_dist{i},animalmodel{animalnum}.lmodel.rep))/2);
            nogovec(animalnum,i) = nanmean((cellfun(@(x) x.part{1}.nogo.lda_norm_dist{i} + x.part{2}.nogo.lda_norm_dist{i},animalmodel{animalnum}.lmodel.rep))/2);
            %             shgovec(animalnum,i) = nanmean((cellfun(@(x) x.part{1}.go.lda_norm_dist{i} + x.part{2}.go.lda_norm_dist{i},sh.animalmodel{animalnum}.lmodel.rep))/2);
            %             shnogovec(animalnum,i) = nanmean((cellfun(@(x) x.part{1}.nogo.lda_norm_dist{i} + x.part{2}.nogo.lda_norm_dist{i},sh.animalmodel{animalnum}.lmodel.rep))/2);
            shgovec(animalnum,i,:,:) = [(cellfun(@(x) x.part{1}.go.lda_norm_dist{i} ,sh.animalmodel{animalnum}.lmodel.rep));(cellfun(@(x) x.part{2}.go.lda_norm_dist{i} ,sh.animalmodel{animalnum}.lmodel.rep))]';
            shnogovec(animalnum,i,:,:) = [(cellfun(@(x) x.part{1}.nogo.lda_norm_dist{i} ,sh.animalmodel{animalnum}.lmodel.rep));(cellfun(@(x) x.part{2}.nogo.lda_norm_dist{i} ,sh.animalmodel{animalnum}.lmodel.rep))]';
            
        elseif strcmp(distMeasure,'mu-norm')
            govec(animalnum,i) = nanmean((cellfun(@(x) x.part{1}.go.mu_dist_norm{i} + x.part{2}.go.mu_dist_norm{i},animalmodel{animalnum}.lmodel.rep))/2);
            nogovec(animalnum,i) = nanmean((cellfun(@(x) x.part{1}.nogo.mu_dist_norm{i} + x.part{2}.nogo.mu_dist_norm{i},animalmodel{animalnum}.lmodel.rep))/2);
            shgovec(animalnum,i,:,:) = [(cellfun(@(x) x.part{1}.go.mu_dist_norm{i} ,sh.animalmodel{animalnum}.lmodel.rep));(cellfun(@(x) x.part{2}.go.mu_dist_norm{i} ,sh.animalmodel{animalnum}.lmodel.rep))]';
            shnogovec(animalnum,i,:,:) = [(cellfun(@(x) x.part{1}.nogo.mu_dist_norm{i} ,sh.animalmodel{animalnum}.lmodel.rep));(cellfun(@(x) x.part{2}.nogo.mu_dist_norm{i} ,sh.animalmodel{animalnum}.lmodel.rep))]';
            
        end
    end
end

shvec = nan(numanimals*200,8);
for i = 1:8
    shvec(:,i) = (reshape(shgovec(:,i,:),1,[])+reshape(shnogovec(:,i,:),1,[]))/2;
end

anova1(govec)
anova1(nogovec)
%%%%%%%%%%%

figure;subplot(2,1,1)
title(distMeasure)
hold on; scatter((1:8),nanmean(govec,1)','.g')
hold on; errorbar((1:8),nanmean(govec,1)',1*nanstd(govec)/sqrt(numanimals),'LineStyle','none','Color','g')
hold on;%bar((1:8)+0.2,nanmean(nogovec,1)','BarWidth',0.4,'FaceColor','r');
hold on;plot(1:8,nanmean(shvec,1),'k-');
if add_chance_ci
hold on;plot(1:8,nanmean(shvec,1)-2*nanstd(shvec),'k--');
hold on;plot(1:8,nanmean(shvec,1)+2*nanstd(shvec),'k--');
end
xlim([0,9])
if strcmp(exptype,'FF')
    ylim([0 4]);
else
    ylim([0 1]);
end

subplot(2,1,2)
hold on; scatter((1:8),nanmean(nogovec,1)','.r')
hold on; errorbar((1:8)+0,nanmean(nogovec,1)',1*nanstd(nogovec)/sqrt(numanimals),'LineStyle','none','Color','r')
hold on;plot(1:8,nanmean(shvec,1),'k-');
if add_chance_ci
hold on;plot(1:8,nanmean(shvec,1)-2*nanstd(shvec),'k--');
hold on;plot(1:8,nanmean(shvec,1)+2*nanstd(shvec),'k--');
end
xlim([0,9])
if strcmp(exptype,'FF')
    ylim([0 4]);
else
    ylim([0 1]);
end
%%%
figure;%bar((1:8)-0.2,nanmean(govec,1)','BarWidth',0.4,'FaceColor','g');
hold on; scatter((1:8)-0.1,nanmean(govec,1)','.g')
hold on; errorbar((1:8)-0.1,nanmean(govec,1)',1*nanstd(govec)/sqrt(numanimals),'LineStyle','none','Color','g')
hold on;%bar((1:8)+0.2,nanmean(nogovec,1)','BarWidth',0.4,'FaceColor','r');
hold on; scatter((1:8)+0.1,nanmean(nogovec,1)','.r')
hold on; errorbar((1:8)+0.1,nanmean(nogovec,1)',1*nanstd(nogovec)/sqrt(numanimals),'LineStyle','none','Color','r')
hold on;plot(1:8,nanmean(shvec,1),'k-');
if add_chance_ci
hold on;plot(1:8,nanmean(shvec,1)-2*nanstd(shvec),'k--');
hold on;plot(1:8,nanmean(shvec,1)+2*nanstd(shvec),'k--');
end
if strcmp(exptype,'FF')
    ylim([0 4]);
else
    ylim([0 1]);
end


pvals = nan(1,8);
for i = 1:8
 pvals(i) = ranksum(govec(:,i),nogovec(:,i));
end
pvals