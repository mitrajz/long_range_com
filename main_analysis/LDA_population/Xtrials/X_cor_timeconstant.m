
%% should make animalmodel files before here, communication maps per animals are loaded from animalmodel (LDA results)

% [in the latest version, no need to load cells, just run (activity is taken from lda model)]
% [in the latest version, estimates of tau are done with jackknife, which is timetaking. 
% For a quick check, comment ci estimation in *_plot_helper.m]

exptype = 'FB';
binsizems = '65';
addshuffle = 0; % default 1, but is unnecessary for most things
nrep = '100';
xstyle = '2'; % default 2
ldaormu = 1; %1: using lda as com vectors 1: using mu as com vectors - activity vectors would be unchanged
peranimal = 0;% if 1: calculates tau per animal as well as the combied one.
saveres = 0; % if 1: automatically makes timestamped folder and saves figures and tau estimates.
normalizeorsource = 2; % if 1: raw and normalized targetactivity, if 2:raw target and raw source
fitmethod = 'T'; % L:levenberg, T:trust region
% no need to load cells. animalcells contain all
% to visualize shuffle cormats, jst add 'shuffle_' to the beginning of the
% filename and set addshuffle to 0:
%ldafilename = ['shuffle_',exptype,'lmodel_25Th_',binsizems,'ms_X_style',xstyle,'_nrep',nrep,'.mat'];

if normalizeorsource == 1
    n1type = 'withnorm';
elseif normalizeorsource == 2
    n1type = 'withsource';
end

preSub = '_preSub2'; % options: '' (default), '_preSub', and '_preSub2'

savdir = '/mnt/data/Mitra/figs/P2_L/bothdircombined/2020_02_24/LDA/version2/cosinesimilaritymats_and_timeconstants/saveres';
ldafilename = [exptype,'lmodel_25Th_',binsizems,'ms_X_style',xstyle,preSub,'_nrep',nrep,'.mat'];
filename = fullfile('/mnt/data/Mitra/figs/P2_L/bothdircombined/2020_02_24/LDA/version2',ldafilename);
if addshuffle
    SH = load(fullfile('/mnt/data/Mitra/figs/P2_L/bothdircombined/2020_02_24/LDA/version2',['shuffle_',ldafilename]));
end
timescalems = 65; % used just for labeling xaxis - 
nlags = 8;
% keep 0, except if you want to combine FF and FB for baseline V1 and LM
% activity. In which case keep 1. And change round from 1 to 2 when running
% second round
combined = 0; %%% should keep 0, each cell type for one experiment type (V1 only from FB exps)
round = 1;


load(filename);

allanimals = find(cellfun(@(x) isfield(x,'lmodel'),animalmodel))
numanimals =length(allanimals);
%%%%% exp fit parameters

if strcmp(fitmethod,'T')
        fo = fitoptions('Method','NonlinearLeastSquares','Lower',[-inf,-inf,0],'Upper',[inf,inf,inf],'StartPoint',[0.5,1,1000],'MaxIter',1000);    
%     fo = fitoptions('Method','NonlinearLeastSquares','Lower',[-inf,-inf,0],'Upper',[inf,inf,inf],'StartPoint',[0.5,1,1000],'MaxIter',1000);    
    fo_lin = fitoptions('Method','NonlinearLeastSquares','Lower',[-inf,-inf],'Upper',[inf,inf],'StartPoint',[0 0 ],'MaxIter',1000);
elseif strcmp(fitmethod,'L')
    fo = fitoptions('Method','NonlinearLeastSquares','StartPoint',[0.5,1,2000],'MaxIter',2000,...
        'Algorithm','Levenberg-Marquardt'); 
    fo_lin = fitoptions('Method','NonlinearLeastSquares','StartPoint',[0,0],'MaxIter',2000,...
        'Algorithm','Levenberg-Marquardt');    
else
    error('no method spcified')
end

ftype = fittype('A*(exp(-x/tau)+B)','independent','x','options',fo);
ftype_lin = fittype('A*x+B','independent','x','options',fo_lin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% targetcell and notargetcell should have been specified before here
% load, clean and remove before here (aggregate analysis)

% plotting params
if strcmp(exptype,'FB')
    celltypename = 'V1';
    notargetcelltypename = 'LM';
elseif strcmp(exptype,'FF')
    celltypename = 'LM';
    notargetcelltypename = 'V1';
end
scAlpha=0.5;
xjitter=0.1;
sc = 0; % whether to overlap scattered points




if combined
    if round == 1
    allanimalsn0cormat_go = [];
    allanimalsn0cormat_nogo = [];
    allanimalsn1cormat_go = [];
    allanimalsn1cormat_nogo = [];
    allanimalscomcormat_go = [];
    allanimalscomcormat_nogo = [];
    go_com_tc = nan(5*numanimals*nlags,nlags);
    go_n0_tc = nan(5*numanimals*nlags,nlags);
    go_n1_tc = nan(5*numanimals*nlags,nlags);
    nogo_com_tc = nan(5*numanimals*nlags,nlags);
    nogo_n0_tc = nan(5*numanimals*nlags,nlags);
    nogo_n1_tc = nan(5*numanimals*nlags,nlags);
    
    go_currentrow = 0;
    nogo_currentrow = 0;
    end
else
    allanimalsn0cormat_go = [];
    allanimalsn0cormat_nogo = [];
    allanimalsn1cormat_go = [];
    allanimalsn1cormat_nogo = [];
    allanimalscomcormat_go = [];
    allanimalscomcormat_nogo = [];
    go_com_tc = nan(numanimals*nlags,nlags);
    SH_go_com_tc = nan(numanimals*nlags,nlags);
    go_n0_tc = nan(numanimals*nlags,nlags);
    go_n1_tc = nan(numanimals*nlags,nlags);
    nogo_com_tc = nan(numanimals*nlags,nlags);
    SH_nogo_com_tc = nan(numanimals*nlags,nlags);
    nogo_n0_tc = nan(numanimals*nlags,nlags);
    nogo_n1_tc = nan(numanimals*nlags,nlags);
    
    a_go_com_tc = nan(nlags,nlags);
    a_go_n0_tc = nan(nlags,nlags);
    a_go_n1_tc = nan(nlags,nlags);
    a_nogo_com_tc = nan(nlags,nlags);
    a_nogo_n0_tc = nan(nlags,nlags);
    a_nogo_n1_tc = nan(nlags,nlags);
    
    go_currentrow = 0;
    nogo_currentrow = 0;
end
animaltau_com_go = {};
animaltau_com_nogo = {};
animaltau_n0_go = {};
animaltau_n0_nogo = {};
animaltau_n1_go = {};
animaltau_n1_nogo = {};

tau_com_go = {};
tau_com_nogo = {};
tau_n0_go = {};
tau_n0_nogo = {};
tau_n1_go = {};
tau_n1_nogo = {};

for animalnum = allanimals
    
    animal_go_currentrow = 0;
    animal_nogo_currentrow = 0;
    
    % go: 
       
    accormat_n0 =animalmodel{animalnum}.lmodel.acmatgo_n0;
    if normalizeorsource == 1
        accormat_n1 = animalmodel{animalnum}.lmodel.acmatgo_n1;
    elseif normalizeorsource == 2
        accormat_n1 = animalmodel{animalnum}.lmodel.acmatgo_notarget_n0;
    end
    if ldaormu == 2
        Comcormat = animalmodel{animalnum}.lmodel.Cmatgo_mu;
    else
        Comcormat = animalmodel{animalnum}.lmodel.Cmatgo_lda;
    end
    
    if addshuffle 
      
        if ldaormu == 2
              SH_Comcormat = SH.animalmodel{animalnum}.lmodel.Cmatgo_mu;        
        else
              SH_Comcormat = SH.animalmodel{animalnum}.lmodel.Cmatgo_lda;        
        end
        
    end
    
    allanimalsn0cormat_go = cat(3,allanimalsn0cormat_go,accormat_n0);
    allanimalsn1cormat_go = cat(3,allanimalsn1cormat_go,accormat_n1);
    allanimalscomcormat_go = cat(3,allanimalscomcormat_go,Comcormat);

    
    for ii=1:nlags
        go_currentrow = go_currentrow + 1;
        animal_go_currentrow = animal_go_currentrow + 1;
        for jj =1:nlags
            if (jj-ii) >= 0
                go_com_tc(go_currentrow,jj-ii+1) = Comcormat(ii,jj);
                go_n0_tc(go_currentrow,jj-ii+1) = accormat_n0(ii,jj);
                go_n1_tc(go_currentrow,jj-ii+1) = accormat_n1(ii,jj);
                if addshuffle
                    SH_go_com_tc(go_currentrow,jj-ii+1) = SH_Comcormat(ii,jj);
                end
                a_go_com_tc(animal_go_currentrow,jj-ii+1) = Comcormat(ii,jj);
                a_go_n0_tc(animal_go_currentrow,jj-ii+1) = accormat_n0(ii,jj);
                a_go_n1_tc(animal_go_currentrow,jj-ii+1) = accormat_n1(ii,jj);
                
            end
        end
    end
    
    
    % nogo: 
    accormat_n0 =animalmodel{animalnum}.lmodel.acmatnogo_n0;
    if normalizeorsource == 1
        accormat_n1 = animalmodel{animalnum}.lmodel.acmatnogo_n1;
    elseif normalizeorsource == 2
         accormat_n1 = animalmodel{animalnum}.lmodel.acmatnogo_notarget_n0;
    end
    if ldaormu == 2
        Comcormat = animalmodel{animalnum}.lmodel.Cmatnogo_mu;
    else
        Comcormat = animalmodel{animalnum}.lmodel.Cmatnogo_lda;
    end
    
    if addshuffle 
        if ldaormu == 2
            SH_Comcormat = SH.animalmodel{animalnum}.lmodel.Cmatnogo_mu;
        else
            SH_Comcormat = SH.animalmodel{animalnum}.lmodel.Cmatnogo_lda;
        end
    end
    
    allanimalsn0cormat_nogo = cat(3,allanimalsn0cormat_nogo,accormat_n0);
    allanimalsn1cormat_nogo = cat(3,allanimalsn1cormat_nogo,accormat_n1);
    allanimalscomcormat_nogo = cat(3,allanimalscomcormat_nogo,Comcormat);

    
     for ii=1:nlags
        nogo_currentrow = nogo_currentrow + 1;
        animal_nogo_currentrow = animal_nogo_currentrow + 1;
        for jj =1:nlags
            if (jj-ii) >= 0
                nogo_com_tc(nogo_currentrow,jj-ii+1) = Comcormat(ii,jj);
                nogo_n0_tc(nogo_currentrow,jj-ii+1) = accormat_n0(ii,jj);
                nogo_n1_tc(nogo_currentrow,jj-ii+1) = accormat_n1(ii,jj);
                 if addshuffle
                    SH_nogo_com_tc(nogo_currentrow,jj-ii+1) = SH_Comcormat(ii,jj);
                 end
                a_nogo_com_tc(animal_nogo_currentrow,jj-ii+1) = Comcormat(ii,jj);
                a_nogo_n0_tc(animal_nogo_currentrow,jj-ii+1) = accormat_n0(ii,jj);
                a_nogo_n1_tc(animal_nogo_currentrow,jj-ii+1) = accormat_n1(ii,jj);
            end
        end
     end
    
    if peranimal
        figure;[animaltau_com_go{end+1},~]=X_cor_timeconstant_plot_helper(a_go_com_tc,'g',scAlpha,xjitter,sc,ftype,fo,timescalems,'exp'); %untrained: (go_com_tc+nogo_com_tc)/2,'k'
        hold on;[animaltau_com_nogo{end+1},~]=X_cor_timeconstant_plot_helper(a_nogo_com_tc,'r',scAlpha,-xjitter,sc,ftype,fo,timescalems,'exp');
        title(['animal',num2str(animalnum),exptype,'-',binsizems,'ms'])
        
        figure;[animaltau_n0_go{end+1},~]=X_cor_timeconstant_plot_helper(a_go_n0_tc,'g',scAlpha,xjitter,sc,ftype_lin,fo_lin,timescalems,'lin');
        hold on;[animaltau_n0_nogo{end+1},~]=X_cor_timeconstant_plot_helper(a_nogo_n0_tc,'r',scAlpha,-xjitter,sc,ftype_lin,fo_lin,timescalems,'lin');
        title(['animal',num2str(animalnum),celltypename,'-raw-',binsizems,'ms'])
          
        % skip this for most purposes. Timetaking 
         figure;[animaltau_n1_go{end+1},~]=X_cor_timeconstant_plot_helper(a_go_n1_tc,'g',scAlpha,xjitter,sc,ftype_lin,fo_lin,timescalems,'lin');
         hold on;[animaltau_n1_nogo{end+1},~]=X_cor_timeconstant_plot_helper(a_nogo_n1_tc,'r',scAlpha,-xjitter,sc,ftype_lin,fo_lin,timescalems,'lin');
         if normalizeorsource == 1
            title(['animal',num2str(animalnum),celltypename,'-normalized-',binsizems,'ms'])
         elseif normalizeorsource == 2
             title(['animal',num2str(animalnum),notargetcelltypename,'-raw-',binsizems,'ms'])
         end
      

    end
     
    
end

%%% plotts:
% 2 styles of plots in helper: 1- data points, mean, and mean ci - you see the difference
%                              2- exponential fits with 0 intercept. shaded
%                              bars are prediction intervals of the
%                              function, fitting all simultaneously

clear gof


figure;[tau1,gof(1)]=X_cor_timeconstant_plot_helper(go_n0_tc,'g',scAlpha,xjitter,sc,ftype_lin,fo_lin,timescalems,'lin');
hold on;[tau2,gof(2)]=X_cor_timeconstant_plot_helper(nogo_n0_tc,'r',scAlpha,-xjitter,sc,ftype_lin,fo_lin,timescalems,'lin');
title([celltypename,'-raw-',binsizems,'ms'])
tau_n0_go{end+1} = tau1;
tau_n0_nogo{end+1} = tau2;

figure;[tau3,gof(3)]=X_cor_timeconstant_plot_helper(go_n1_tc,'g',scAlpha,xjitter,sc,ftype_lin,fo_lin,timescalems,'lin');
hold on;[tau4,gof(4)]=X_cor_timeconstant_plot_helper(nogo_n1_tc,'r',scAlpha,-xjitter,sc,ftype_lin,fo_lin,timescalems,'lin');
if normalizeorsource == 1
    title([celltypename,'-normalized-',binsizems,'ms'])
elseif normalizeorsource == 2
    title([notargetcelltypename,'-raw-',binsizems,'ms'])
end
tau_n1_go{end+1} = tau3;
tau_n1_nogo{end+1} = tau4;

figure;[tau5,gof(5)]=X_cor_timeconstant_plot_helper(go_com_tc,'g',scAlpha,xjitter,sc,ftype,fo,timescalems,'exp'); %untrained: (go_com_tc+nogo_com_tc)/2,'k'
hold on;[tau6,gof(6)]=X_cor_timeconstant_plot_helper(nogo_com_tc,'r',scAlpha,-xjitter,sc,ftype,fo,timescalems,'exp');
tau_com_go{end+1} = tau5;
tau_com_nogo{end+1} = tau6;
title([exptype,'-',binsizems,'ms'])
if addshuffle
hold on;[~,~]=X_cor_timeconstant_plot_helper((SH_go_com_tc+SH_nogo_com_tc)/2,'k',scAlpha,xjitter,sc,ftype,fo,timescalems,'exp'); %untrained: (go_com_tc+nogo_com_tc)/2,'k'
%hold on;[~,~]=X_cor_timeconstant_plot_helper(SH_nogo_com_tc,'k',scAlpha,-xjitter,sc,ftype,timescalems);
end
%%%% com exp fit stats
[taugo,taunogo,taudif] = expstats(go_com_tc,nogo_com_tc,ftype,timescalems);


% plot tau
figure;
bar([1 2 3]-0.1,[tau1.val tau3.val tau5.val],0.2,'FaceColor','g','EdgeColor','none','FaceAlpha',0.5)
hold on
errorbar([1 2 3]-0.1,[tau1.val tau3.val tau5.val],...
    [tau1.ci tau3.ci tau5.ci],'.','Color','g')
hold on
bar([1 2 3]+0.1,[tau2.val tau4.val tau6.val],0.2,'FaceColor','r','EdgeColor','none','FaceAlpha',0.5)
errorbar([1 2 3]+0.1,[tau2.val tau4.val tau6.val],...
    [tau2.ci tau4.ci tau6.ci],'.','Color','r')
xlim([0 4])
set(gca,'XTick',1:3)
if normalizeorsource == 1
    set(gca,'XTickLabel',{[celltypename,'raw'],[celltypename,'norm'],'com'})
elseif normalizeorsource == 2
    set(gca,'XTickLabel',{[celltypename,'raw'],[notargetcelltypename,'raw'],'com'})
end
title('tau-ci')

%%%%% V1/LM cormats

figure;
s1=subplot(1,2,1);a1=imagesc(nanmean(allanimalsn0cormat_go,3));
s2=subplot(1,2,2);a2=imagesc(nanmean(allanimalsn0cormat_nogo,3));
s1.Title.String = (sprintf('go-raw-N = %d animals',size(allanimalsn0cormat_go,3)));
s1.Colormap = colormap(jet(64));% redblue
set(s1,'CLim',[-1 1])
colorbar(a1.Parent)
s2.Title.String = (sprintf('nogo-raw-N = %d animals',size(allanimalsn0cormat_nogo,3)));
s2.Colormap = colormap(jet(64));% redblue
set(s2,'CLim',[-1 1])
colorbar(a2.Parent)


figure;
s1=subplot(1,2,1);a1=imagesc(nanmean(allanimalsn1cormat_go,3));
s2=subplot(1,2,2);a2=imagesc(nanmean(allanimalsn1cormat_nogo,3));
if normalizeorsource == 1
    s1.Title.String = (sprintf('go-norm-N = %d animals',size(allanimalsn1cormat_go,3)));
elseif normalizeorsource == 2
    s1.Title.String = (sprintf('go-raw-source-N = %d animals',size(allanimalsn1cormat_go,3)));
end
s1.Colormap = colormap(jet(64));% redblue
set(s1,'CLim',[-1 1])
colorbar(a1.Parent)
if normalizeorsource == 1
    s2.Title.String = (sprintf('nogo-norm-N = %d animals',size(allanimalsn1cormat_nogo,3)));
elseif normalizeorsource == 2
    s2.Title.String = (sprintf('nogo-raw-source-N = %d animals',size(allanimalsn1cormat_nogo,3)));
end
s2.Colormap = colormap(jet(64));% redblue
set(s2,'CLim',[-1 1])
colorbar(a2.Parent)

%%%%% com cormats

figure;
s1=subplot(1,2,1);a1=imagesc(nanmean(allanimalscomcormat_go,3));
s2=subplot(1,2,2);a2=imagesc(nanmean(allanimalscomcormat_nogo,3));
s1.Title.String = (sprintf('go-%s-N = %d animals',exptype,size(allanimalscomcormat_go,3)));
s1.Colormap = colormap(jet(64));% redblue
set(s1,'CLim',[-0.7 0.7])
colorbar(a1.Parent)
s2.Title.String = (sprintf('nogo-%s-N = %d animals',exptype,size(allanimalscomcormat_go,3)));
s2.Colormap = colormap(jet(64));% redblue
set(s2,'CLim',[-0.7 0.7])
colorbar(a2.Parent)
%%%% saving

if saveres
    cd(savdir)
    foldername = [exptype,binsizems,'Xstyle',xstyle,'nrep',nrep,'_',n1type,'_','ldaormu',num2str(ldaormu),preSub,...
        '_',datestr(now,'yyyy_mm_dd_HH_MM_SS')];
    mkdir(foldername)
    cd(foldername)
    % figures
    allfigs = findobj('Type','figure');
    for i = 1:length(allfigs)
        % titles of all axes
        alltitles = arrayfun(@(x) x.Title.String, allfigs(i).Children,'UniformOutput',0);
        % choose one non-empty title
        if ~isempty(alltitles)
            figname = alltitles{find(cellfun(@(x) numel(x), alltitles),1)};
        else
            figname = 'unnamed';
        end
        saveas(allfigs(i),[figname,'.fig']);        
    end
    % variables
    save('alltaus.mat',...
        'animaltau_com_go','animaltau_com_nogo','animaltau_n0_go','animaltau_n0_nogo','animaltau_n1_go','animaltau_n1_nogo',...
        'tau_com_go','tau_com_nogo','tau_n0_go','tau_n0_nogo','tau_n1_go','tau_n1_nogo')
end

