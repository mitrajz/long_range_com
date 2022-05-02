%%%%%%%% consistency plots 
exptype = 'FF';
binsizems = '150';
Xstyle = 1; 
nrep = 100; %


filename = [exptype,'lmodel_25Th_',binsizems,'ms_X_style',num2str(Xstyle),'_nrep',num2str(nrep),'.mat'];
D = load(fullfile('/mnt/data/Mitra/figs/P2_L/bothdircombined/2020_02_24/LDA/version2',filename));
SH = load(fullfile('/mnt/data/Mitra/figs/P2_L/bothdircombined/2020_02_24/LDA/version2',['shuffle_',filename]));
includedanimals = find(cellfun(@(x) isfield(x,'lmodel'),D.animalmodel))


%%%%%%%%
[Cmatgo_lda,Cmatnogo_lda,Crepgo_lda,Crepnogo_lda,C_go_lda,C_nogo_lda] = X_prepareCmatrices_lda(D.animalmodel,includedanimals);
[Cmatgo_mu,Cmatnogo_mu,Crepgo_mu,Crepnogo_mu,C_go_mu,C_nogo_mu] = X_prepareCmatrices_mu(D.animalmodel,includedanimals);

[~,~,~,~,SH_C_go_lda,SH_C_nogo_lda] = X_prepareCmatrices_lda(SH.animalmodel,includedanimals);
[~,~,~,~,SH_C_go_mu,SH_C_nogo_mu] = X_prepareCmatrices_mu(SH.animalmodel,includedanimals);
SH_C_lda = (SH_C_go_lda+SH_C_nogo_lda)/2;
SH_C_mu = (SH_C_go_mu+SH_C_nogo_mu)/2;

figure
%%%%%%% histograms

s=subplot(2,2,1); 
histogram(C_go_lda(:,1),-0.2:0.1:1,'FaceColor','g','Normalization','pdf','EdgeAlpha',0);hold on;
histogram(C_nogo_lda(:,1),-0.2:0.1:1,'FaceColor','r','Normalization','pdf','EdgeAlpha',0);hold on;
histogram(SH_C_lda(:,1),-0.2:0.1:1,'FaceColor','k','Normalization','pdf','EdgeAlpha',0);hold on;
s.Title.String = sprintf('lda-consistency-Nrep=%d,Xstyle=%d',D.animalmodel{1}.lmodel.nrep,D.animalmodel{1}.lmodel.Xstyle);

s=subplot(2,2,3); 
histogram(C_go_mu(:,1),-0.2:0.1:1,'FaceColor','g','Normalization','pdf','EdgeAlpha',0);hold on;
histogram(C_nogo_mu(:,1),-0.2:0.1:1,'FaceColor','r','Normalization','pdf','EdgeAlpha',0);hold on;
histogram(SH_C_mu(:,1),-0.2:0.1:1,'FaceColor','k','Normalization','pdf','EdgeAlpha',0);hold on;
s.Title.String = sprintf('mu-consistency-Nrep=%d,Xstyle=%d',D.animalmodel{1}.lmodel.nrep,D.animalmodel{1}.lmodel.Xstyle);


%%%%%% errorbars: they are ste and not std
s=subplot(2,2,2);
errorbar(0.9,nanmean(C_go_lda(:,1),1),2*nanstd(C_go_lda(:,1))/sqrt(numel(find(~isnan(C_go_lda(:,1))))),'g')
hold on;
errorbar(1.1,nanmean(C_nogo_lda(:,1),1),2*nanstd(C_nogo_lda(:,1))/sqrt(numel(find(~isnan(C_nogo_lda(:,1))))),'r');
hold on;
scatter(0.9,nanmean(C_go_lda(:,1),1),'.g')
hold on;
scatter(1.1,nanmean(C_nogo_lda(:,1),1),'.r');

hold on
errorbar(1.3,nanmean(SH_C_lda(:,1),1),2*nanstd(SH_C_lda(:,1))/sqrt(numel(find(~isnan(SH_C_lda(:,1))))),'k')
hold on;
scatter(1.3,nanmean(SH_C_lda(:,1),1),'.k')


text(0.5,0.8,sprintf('ranksum p = %d',ranksum(C_go_lda(:,1),C_nogo_lda(:,1))))
xlim([0.5 1.5])
s.Title.String = sprintf('lda-consistency-Nrep=%d,Xstyle=%d',D.animalmodel{1}.lmodel.nrep,D.animalmodel{1}.lmodel.Xstyle);
ylim([-0.5 1])

s=subplot(2,2,4);
errorbar(0.9,nanmean(C_go_mu(:,1),1),2*nanstd(C_go_mu(:,1))/sqrt(numel(find(~isnan(C_go_mu(:,1))))),'g')
hold on;
errorbar(1.1,nanmean(C_nogo_mu(:,1),1),2*nanstd(C_nogo_mu(:,1))/sqrt(numel(find(~isnan(C_nogo_mu(:,1))))),'r');
hold on;
scatter(0.9,nanmean(C_go_mu(:,1),1),'.g')
hold on;
scatter(1.1,nanmean(C_nogo_mu(:,1),1),'.r');

hold on
errorbar(1.3,nanmean(SH_C_mu(:,1),1),2*nanstd(SH_C_mu(:,1))/sqrt(numel(find(~isnan(SH_C_mu(:,1))))),'k')
hold on;
scatter(1.3,nanmean(SH_C_mu(:,1),1),'.k')

xlim([0.5 1.5])
s.Title.String = sprintf('mu-consistency-Nrep=%d,Xstyle=%d',D.animalmodel{1}.lmodel.nrep,D.animalmodel{1}.lmodel.Xstyle);
text(0.5,0.8,sprintf('ranksum p = %d',ranksum(C_go_mu(:,1),C_nogo_mu(:,1))))
ylim([-0.5 1])


