% loads lda, and calculates the angle between [lda(t)] and [activity(t+1) -
% activity(t)]. 
%
exptype = 'FB';
cd('/mnt/data/Mitra/figs/P2_L/bothdircombined/2020_02_24/LDA/version1')
load([exptype,'lmodel_25Th_shufflebsbs.mat'])
load(['cells_',exptype,'_pl20_an150_lw20_exG1_onlyC1_onlyS0_plstyle1.mat']);
cleanupcriteria = struct;
cleanupcriteria.doclean = 1; % 1 or 0, perform any cleaning or not. If zero does nothing
cleanupcriteria.significance = nan; % significance of activity instead of absolute value
cleanupcriteria.smb_activity = 2.5; % in Hz

%%%%

cellremovecriteria.doclean = 1; % 1 or 0, perform any cleaning or not. If zero does nothing
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

%%% clean and remove cells

%load(cellfiletoload);

[LMcells,V1cells,Behcells,params] = cleanupcells(cleanupcriteria,LMcells,V1cells,Behcells,params);
[cell_keep_ind,LMcells,V1cells,Behcells,params] = removecells(cell_keep_ind,cellremovecriteria,LMcells,V1cells,Behcells,params);

%%%%%%%%%%%%%%
nrep = 100;
if strcmp(exptype,'FF')
    targetcell = LMcells;
    nanimals = 7;
elseif strcmp(exptype,'FB')
    targetcell = V1cells;
    nanimals = 6;
end
nlags = 8;
angs_go = nan(nanimals,7); %nanimals*nlags-1
angs_nogo = nan(nanimals,7);
angs_go_sh = nan(nanimals*nrep,7);
angs_nogo_sh  = nan(nanimals*nrep,7);
ldaadd = 1; % 0 or 1

for animalnum = 1:nanimals
    
    animalcellind=find(cellfun(@(x) x.simulcode,targetcell) == animalnum );
    
    
    randeq = 0; % if zero, takes last. Seed is fixed  for random sampling
    targetcell(animalcellind) = equalizegonogotrials(targetcell(animalcellind),0,0,randeq);% frnormalize
    
    for lag =1:(nlags-1)
        
        meandifnogo = nanmean(cell2mat(cellfun(@(x) x.laAbs.nogo(lag+1),targetcell(animalcellind))),1) - ...
            nanmean(cell2mat(cellfun(@(x) x.laAbs.nogo(lag),targetcell(animalcellind))),1);
        
        meandifgo = nanmean(cell2mat(cellfun(@(x) x.laAbs.go(lag+1),targetcell(animalcellind))),1) - ...
            nanmean(cell2mat(cellfun(@(x) x.laAbs.go(lag),targetcell(animalcellind))),1);
        

        meandifgo(find(isnan(meandifgo))) = 0;
        meandifnogo(find(isnan(meandifnogo))) = 0;
        
        
        LDAvecgo = animalmodel{animalnum}.lmodel.wt(lag+ldaadd,:);
        LDAvecnogo = animalmodel{animalnum}.lmodel.wt(lag+8+ldaadd,:);
        
        angs_go(animalnum,lag) = ((meandifgo./norm(meandifgo) )* LDAvecgo');
        angs_nogo(animalnum,lag) = ((meandifnogo./norm(meandifnogo) )* LDAvecnogo');
        
        % now add for shuffle:
        for rep = 1:nrep
            % bs population go at  time lag
            population_go_lag =  ...
                animalmodel{animalnum}.lmodel.shuffle_go{lag}.X{rep}(find(animalmodel{animalnum}.lmodel.shuffle_go{lag}.Y{rep}==0),:);
            
            population_go_lagp1 =  ...
                animalmodel{animalnum}.lmodel.shuffle_go{lag+1}.X{rep}(find(animalmodel{animalnum}.lmodel.shuffle_go{lag+1}.Y{rep}==0),:);
            
            meandifgo_sh = nanmean(population_go_lagp1,1) - nanmean(population_go_lag,1);
            
            population_nogo_lag =  ...
                animalmodel{animalnum}.lmodel.shuffle_nogo{lag}.X{rep}(find(animalmodel{animalnum}.lmodel.shuffle_nogo{lag}.Y{rep}==0),:);
            
            population_nogo_lagp1 =  ...
                animalmodel{animalnum}.lmodel.shuffle_nogo{lag+1}.X{rep}(find(animalmodel{animalnum}.lmodel.shuffle_nogo{lag+1}.Y{rep}==0),:);
            
            meandifnogo_sh = nanmean(population_nogo_lagp1,1) - nanmean(population_nogo_lag,1);
            

            meandifgo_sh(find(isnan(meandifgo_sh))) = 0;
            meandifnogo_sh(find(isnan(meandifnogo_sh))) = 0;
            
            angs_go_sh((animalnum-1)*nrep+rep,lag) = ((meandifgo_sh./norm(meandifgo_sh) )* animalmodel{animalnum}.lmodel.shuffle_go{lag+ldaadd}.vec_n(:,rep));
            angs_nogo_sh((animalnum-1)*nrep+rep,lag) = ((meandifnogo_sh./norm(meandifnogo_sh) ) * animalmodel{animalnum}.lmodel.shuffle_nogo{lag+ldaadd}.vec_n(:,rep));
            
        end
        
        
    end
end

figure
histogram(abs(reshape(angs_go,1,[])),-1:0.1:1,'faceColor','g','EdgeAlpha',0,'FaceAlpha',0.2,'Normalization','pdf');
hold on
histogram(abs(reshape(angs_nogo,1,[])),-1:0.1:1,'faceColor','r','EdgeAlpha',0,'FaceAlpha',0.2,'Normalization','pdf');
hold on;
histogram(abs(reshape([angs_go_sh,angs_nogo_sh],1,[])),-1:0.1:1,'faceColor',[0.5 0.5 0.5],'EdgeAlpha',0,'FaceAlpha',0.2,'Normalization','pdf');
hold on;
histogram(abs(reshape(angs_go,1,[])),-1:0.1:1,'EdgeColor','g','FaceAlpha',0,'Normalization','pdf','DisplayStyle','stairs');
hold on
histogram(abs(reshape(angs_nogo,1,[])),-1:0.1:1,'EdgeColor','r','FaceAlpha',0,'Normalization','pdf','DisplayStyle','stairs');
hold on;
histogram(abs(reshape([angs_go_sh,angs_nogo_sh],1,[])),-1:0.1:1,'EdgeColor',[0.5 0.5 0.5],'FaceAlpha',0,'Normalization','pdf','DisplayStyle','stairs');

pgovsshuffle = ranksum(abs(reshape(angs_go,1,[])),abs(reshape([angs_go_sh,angs_nogo_sh],1,[])));
pnogovsshuffle = ranksum(abs(reshape(angs_nogo,1,[])),abs(reshape([angs_go_sh,angs_nogo_sh],1,[])));

text(-1,1,sprintf('go vs shuffle p value\n%d',pgovsshuffle));
text(-1,0.5,sprintf('nogo vs shuffle p value\n%d',pnogovsshuffle));


% errorbars are 1*sem and std for shuffle
figure;
errorbar(nanmean(angs_go,1),nanstd(angs_go)/sqrt(size(angs_go,1)),'g');
hold on;
errorbar(nanmean(angs_nogo,1),nanstd(angs_nogo)/sqrt(size(angs_nogo,1)),'r');
hold on;
errorbar(nanmean(angs_go_sh,1),nanstd(angs_go_sh),'k');
hold on;
errorbar(nanmean(angs_nogo_sh,1),nanstd(angs_nogo_sh),'k');
