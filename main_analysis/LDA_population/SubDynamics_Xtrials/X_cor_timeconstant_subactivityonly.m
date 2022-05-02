% simplified version of tc, only for activity, fast, slow, and all the
% cells
% importantly ff and fb experiments are combined here
% 65 ms
% first file source, second file target. 
% ff and fb exps combined
% change ratio to select all cells, slow, or fast cells
% change the order: ff,fb or fb,ff to change V1 or LM
% first, run X_lda_new_subdynamics.m, files are saved with a dummy pefix.
% These files are then loaded here.
%% should make animalmodel files before here, communication maps per animals are loaded from animalmodel (LDA results)
% to run, set exp type, load cells file, aggregate analysis - clean and remove, then run this section 
% Here, normalize is actually the normalized activity (and not saved per
% animal). It is not like X_cor_timeconstant, in which it can be the source
% area. This script produces activity dynamics for 1 area at a time (V1 or
% LM), pooling together both ff and fb expriments. if ratio nan, all cells,
% other wise a proportion of fastest or slowest cells.

binsizems = '65';
nrep = '100';
xstyle = '2';
ratio = nan;%'0.1'; %nan: all cells or '0.2' or '0.8'
exporder ={'FB','FF'};
fitmethod = 'T'; % L:levenberg, T:trust region
peranimal = 1;% 0 or 1 - 1 for slope statistics
saveres = 1; % 0 or 1 - 1 for slope statistics

savdir = '/mnt/data/Mitra/figs/P2_L/bothdircombined/2020_02_24/LDA/version2/cosinesimilaritymats_and_timeconstants/saveres';
% no need to load cells. animalcells contain all

if isnan(ratio)
    first_ldafilename = [exporder{1},'lmodel_25Th_',binsizems,'ms_X_style',xstyle,'_nrep',nrep,'.mat'];
else
    first_ldafilename = ['dummyforactivity',ratio,exporder{1},'lmodel_25Th_',binsizems,'ms_X_style',xstyle,'_nrep',nrep,'.mat'];
end
first_filename = fullfile('/mnt/data/Mitra/figs/P2_L/bothdircombined/2020_02_24/LDA/version2',first_ldafilename);

if isnan(ratio)
    second_ldafilename = [exporder{2},'lmodel_25Th_',binsizems,'ms_X_style',xstyle,'_nrep',nrep,'.mat'];
else
    second_ldafilename = ['dummyforactivity',ratio,exporder{2},'lmodel_25Th_',binsizems,'ms_X_style',xstyle,'_nrep',nrep,'.mat'];
end
second_filename = fullfile('/mnt/data/Mitra/figs/P2_L/bothdircombined/2020_02_24/LDA/version2',second_ldafilename);

timescalems = 65; % used just for labeling xaxis - 
nlags = 8;

%%%%% exp fit parameters

if strcmp(fitmethod,'T')
    fo = fitoptions('Method','NonlinearLeastSquares','Lower',[-inf,-inf,0],'Upper',[inf,inf,inf],'StartPoint',[0.5,1,1000],'MaxIter',1000);    
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



clear gof
% plotting params
scAlpha=0.5;
xjitter=0.1;
sc = 0; % whether to overlap scattered points

if numel(strfind(first_ldafilename,'FF'))
    celltypename = 'V1';
elseif numel(strfind(first_ldafilename,'FB'))
    celltypename = 'LM';
end

if ~isnan(ratio)
if str2num(ratio)<0.5
    celltypename=[celltypename,'-fast',ratio];
else
    celltypename=[celltypename,'-slow',ratio];
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% targetcell and notargetcell should have been specified before here
% load, clean and remove before here (aggregate analysis)

allanimalsn0cormat_go = [];
allanimalsn0cormat_nogo = [];
allanimalsn1cormat_go = [];
allanimalsn1cormat_nogo = [];

maxnumanimals = 14;
go_n0_tc = nan(maxnumanimals*nlags,nlags);
go_n1_tc = nan(maxnumanimals*nlags,nlags);
nogo_n0_tc = nan(maxnumanimals*nlags,nlags);
nogo_n1_tc = nan(maxnumanimals*nlags,nlags);

a_go_n0_tc = nan(nlags,nlags);
a_go_n1_tc = nan(nlags,nlags);
a_nogo_n0_tc = nan(nlags,nlags);
a_nogo_n1_tc = nan(nlags,nlags);

animaltau_n0_go = {};
animaltau_n0_nogo = {};
animaltau_n1_go = {};
animaltau_n1_nogo = {};


go_currentrow = 0;
nogo_currentrow = 0;

% first file source, second target
load(first_filename);
allanimals = find(cellfun(@(x) isfield(x,'lmodel'),animalmodel))
numanimals =length(allanimals);
for animalnum = allanimals   
    
    animal_go_currentrow = 0;
    animal_nogo_currentrow = 0;
    
    % go:        
    accormat_n0 =animalmodel{animalnum}.lmodel.acmatgo_notarget_n0;
    accormat_n1 = animalmodel{animalnum}.lmodel.acmatgo_notarget_n1;
    allanimalsn0cormat_go = cat(3,allanimalsn0cormat_go,accormat_n0);
    allanimalsn1cormat_go = cat(3,allanimalsn1cormat_go,accormat_n1);    
    for ii=1:nlags
        go_currentrow = go_currentrow + 1;
        animal_go_currentrow = animal_go_currentrow + 1;
        for jj =1:nlags
            if (jj-ii) >= 0         
                go_n0_tc(go_currentrow,jj-ii+1) = accormat_n0(ii,jj);
                go_n1_tc(go_currentrow,jj-ii+1) = accormat_n1(ii,jj);    
                
                a_go_n0_tc(animal_go_currentrow,jj-ii+1) = accormat_n0(ii,jj);
                a_go_n1_tc(animal_go_currentrow,jj-ii+1) = accormat_n1(ii,jj);
            end
        end
    end    
    % nogo: 
    accormat_n0 =animalmodel{animalnum}.lmodel.acmatnogo_notarget_n0;
    accormat_n1 = animalmodel{animalnum}.lmodel.acmatnogo_notarget_n1;
    allanimalsn0cormat_nogo = cat(3,allanimalsn0cormat_nogo,accormat_n0);
    allanimalsn1cormat_nogo = cat(3,allanimalsn1cormat_nogo,accormat_n1);      
     for ii=1:nlags
        nogo_currentrow = nogo_currentrow + 1;         
        animal_nogo_currentrow = animal_nogo_currentrow + 1;
        for jj =1:nlags
            if (jj-ii) >= 0
                nogo_n0_tc(nogo_currentrow,jj-ii+1) = accormat_n0(ii,jj);
                nogo_n1_tc(nogo_currentrow,jj-ii+1) = accormat_n1(ii,jj);    
                
                a_nogo_n0_tc(animal_nogo_currentrow,jj-ii+1) = accormat_n0(ii,jj);
                a_nogo_n1_tc(animal_nogo_currentrow,jj-ii+1) = accormat_n1(ii,jj);
            end
        end
     end   
     
      if peranimal
         figure;[animaltau_n0_go{end+1},~]=X_cor_timeconstant_plot_helper(a_go_n0_tc,'g',scAlpha,xjitter,sc,ftype_lin,fo_lin,timescalems,'lin');
        hold on;[animaltau_n0_nogo{end+1},~]=X_cor_timeconstant_plot_helper(a_nogo_n0_tc,'r',scAlpha,-xjitter,sc,ftype_lin,fo_lin,timescalems,'lin');
        title(['animal',num2str(animalnum),celltypename,'-raw-',binsizems,'ms'])
          
%         % skip this for most purposes. Timetaking 
%          figure;[animaltau_n1_go{end+1},~]=X_cor_timeconstant_plot_helper(a_go_n1_tc,'g',scAlpha,xjitter,sc,ftype_lin,fo_lin,timescalems,'lin');
%          hold on;[animaltau_n1_nogo{end+1},~]=X_cor_timeconstant_plot_helper(a_nogo_n1_tc,'r',scAlpha,-xjitter,sc,ftype_lin,fo_lin,timescalems,'lin');
%          if normalizeorsource == 1
%             title(['animal',num2str(animalnum),celltypename,'-normalized-',binsizems,'ms'])
%          elseif normalizeorsource == 2
%              title(['animal',num2str(animalnum),notargetcelltypename,'-raw-',binsizems,'ms'])
%          end
    end
     
     
end
% second
load(second_filename);
allanimals = find(cellfun(@(x) isfield(x,'lmodel'),animalmodel))
numanimals =length(allanimals);
for animalnum = allanimals   
    
    animal_go_currentrow = 0;
    animal_nogo_currentrow = 0;
    
    % go:        
    accormat_n0 =animalmodel{animalnum}.lmodel.acmatgo_n0;
    accormat_n1 = animalmodel{animalnum}.lmodel.acmatgo_n1;
    allanimalsn0cormat_go = cat(3,allanimalsn0cormat_go,accormat_n0);
    allanimalsn1cormat_go = cat(3,allanimalsn1cormat_go,accormat_n1);    
    for ii=1:nlags
        go_currentrow = go_currentrow + 1;
        animal_go_currentrow = animal_go_currentrow + 1;
        for jj =1:nlags
            if (jj-ii) >= 0         
                go_n0_tc(go_currentrow,jj-ii+1) = accormat_n0(ii,jj);
                go_n1_tc(go_currentrow,jj-ii+1) = accormat_n1(ii,jj); 
                
                a_go_n0_tc(animal_go_currentrow,jj-ii+1) = accormat_n0(ii,jj);
                a_go_n1_tc(animal_go_currentrow,jj-ii+1) = accormat_n1(ii,jj);
            end
        end
    end    
    % nogo: 
    accormat_n0 =animalmodel{animalnum}.lmodel.acmatnogo_n0;
    accormat_n1 = animalmodel{animalnum}.lmodel.acmatnogo_n1;
    allanimalsn0cormat_nogo = cat(3,allanimalsn0cormat_nogo,accormat_n0);
    allanimalsn1cormat_nogo = cat(3,allanimalsn1cormat_nogo,accormat_n1);      
     for ii=1:nlags
        nogo_currentrow = nogo_currentrow + 1;
        animal_nogo_currentrow = animal_nogo_currentrow + 1;
        for jj =1:nlags
            if (jj-ii) >= 0
                nogo_n0_tc(nogo_currentrow,jj-ii+1) = accormat_n0(ii,jj);
                nogo_n1_tc(nogo_currentrow,jj-ii+1) = accormat_n1(ii,jj);   
                
                 a_nogo_n0_tc(animal_nogo_currentrow,jj-ii+1) = accormat_n0(ii,jj);
                a_nogo_n1_tc(animal_nogo_currentrow,jj-ii+1) = accormat_n1(ii,jj);
            end
        end
     end   
     
      if peranimal
         figure;[animaltau_n0_go{end+1},~]=X_cor_timeconstant_plot_helper(a_go_n0_tc,'g',scAlpha,xjitter,sc,ftype_lin,fo_lin,timescalems,'lin');
        hold on;[animaltau_n0_nogo{end+1},~]=X_cor_timeconstant_plot_helper(a_nogo_n0_tc,'r',scAlpha,-xjitter,sc,ftype_lin,fo_lin,timescalems,'lin');
        title(['animal',num2str(animalnum),celltypename,'-raw-',binsizems,'ms'])
          
%         % skip this for most purposes. Timetaking 
%          figure;[animaltau_n1_go{end+1},~]=X_cor_timeconstant_plot_helper(a_go_n1_tc,'g',scAlpha,xjitter,sc,ftype_lin,fo_lin,timescalems,'lin');
%          hold on;[animaltau_n1_nogo{end+1},~]=X_cor_timeconstant_plot_helper(a_nogo_n1_tc,'r',scAlpha,-xjitter,sc,ftype_lin,fo_lin,timescalems,'lin');
%          if normalizeorsource == 1
%             title(['animal',num2str(animalnum),celltypename,'-normalized-',binsizems,'ms'])
%          elseif normalizeorsource == 2
%              title(['animal',num2str(animalnum),notargetcelltypename,'-raw-',binsizems,'ms'])
%          end
    end
end





%%% plotts:
% 2 styles of plots in helper: 1- data points, mean, and mean ci - you see the difference
%                              2- exponential fits with 0 intercept. shaded
%                              bars are prediction intervals of the
%                              function, fitting all simultaneously

    
    

figure;[tau1,gof(1)]=X_cor_timeconstant_plot_helper(go_n0_tc,'g',scAlpha,xjitter,sc,ftype_lin,fo_lin,timescalems,'lin');
hold on;[tau2,gof(2)]=X_cor_timeconstant_plot_helper(nogo_n0_tc,'r',scAlpha,-xjitter,sc,ftype_lin,fo_lin,timescalems,'lin');
title([celltypename,'-raw-',binsizems,'ms'])

figure;[tau_c,gof(1)]=X_cor_timeconstant_plot_helper((go_n0_tc+nogo_n0_tc)/2,'k',scAlpha,xjitter,sc,ftype_lin,fo_lin,timescalems,'lin');
title([celltypename,'-raw-',binsizems,'ms','-combiend'])


figure;[tau3,gof(3)]=X_cor_timeconstant_plot_helper(go_n1_tc,'g',scAlpha,xjitter,sc,ftype,fo,timescalems,'exp');
hold on;[tau4,gof(4)]=X_cor_timeconstant_plot_helper(nogo_n1_tc,'r',scAlpha,-xjitter,sc,ftype,fo,timescalems,'exp');
title([celltypename,'-normalized-',binsizems,'ms'])


% plot tau
figure;
bar([1 2 ]-0.1,[tau1.val tau3.val ],0.2,'FaceColor','g','EdgeColor','none','FaceAlpha',0.5)
hold on
errorbar([1 2 ]-0.1,[tau1.val tau3.val ],...
    [tau1.ci tau3.ci ],'.','Color','g')
hold on
bar([1 2 ]+0.1,[tau2.val tau4.val],0.2,'FaceColor','r','EdgeColor','none','FaceAlpha',0.5)
errorbar([1 2]+0.1,[tau2.val tau4.val],...
    [tau2.ci tau4.ci],'.','Color','r')
xlim([0 4])
set(gca,'XTick',1:2)
set(gca,'XTickLabel',{[celltypename,'-raw-',binsizems,'ms'],...
    [celltypename,'-normalized-',binsizems,'ms']})

%%%%% V1/LM cormats

figure;
s1=subplot(1,2,1);a1=imagesc(nanmean(allanimalsn0cormat_go,3));
s2=subplot(1,2,2);a2=imagesc(nanmean(allanimalsn0cormat_nogo,3));
s1.Title.String = (sprintf('go,raw,N = %d animals',size(allanimalsn0cormat_go,3)));
s1.Colormap = colormap(jet(64));% redblue
set(s1,'CLim',[-1 1])
colorbar(a1.Parent)
s2.Title.String = (sprintf('nogo,raw,N = %d animals',size(allanimalsn0cormat_nogo,3)));
s2.Colormap = colormap(jet(64));% redblue
set(s2,'CLim',[-1 1])
colorbar(a2.Parent)


figure;
s1=subplot(1,2,1);a1=imagesc(nanmean(allanimalsn1cormat_go,3));
s2=subplot(1,2,2);a2=imagesc(nanmean(allanimalsn1cormat_nogo,3));
s1.Title.String = (sprintf('go,norm,N = %d animals',size(allanimalsn1cormat_go,3)));
s1.Colormap = colormap(jet(64));% redblue
set(s1,'CLim',[-1 1])
colorbar(a1.Parent)
s2.Title.String = (sprintf('nogo,norm,N = %d animals',size(allanimalsn1cormat_nogo,3)));
s2.Colormap = colormap(jet(64));% redblue
set(s2,'CLim',[-1 1])
colorbar(a2.Parent)



if saveres
    cd(savdir)
    foldername = [celltypename,'_Activity_','ff_fb_combined',binsizems,'Xstyle',xstyle,'nrep',nrep,'_',datestr(now,'yyyy_mm_dd_HH_MM_SS')];
    mkdir(foldername)
    cd(foldername)
    % figures
    allfigs = findobj('Type','figure');
    for i = 1:length(allfigs)
        % titles of all axes
        alltitles = arrayfun(@(x) x.Title.String, allfigs(i).Children,'UniformOutput',0);
        % choose one non-empty title
        try
            figname = alltitles{find(cellfun(@(x) numel(x), alltitles),1)};
        catch
            figname = 'unnamed';
        end
        saveas(allfigs(i),[figname,'.fig']);        
    end
    % variables
    save('alltaus.mat',...
        'animaltau_n0_go','animaltau_n0_nogo')
end

%% saving figures
% while numel(findobj('Type','Figure'))
%     f=gcf;
%     f.Children.Title.String;
%     saveas(f,[f.Children.Title.String,'.fig'])    
%     close(f);
% end
