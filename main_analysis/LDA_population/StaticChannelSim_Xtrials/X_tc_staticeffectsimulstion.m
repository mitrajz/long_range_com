% time constants for the control animal mdoels
% load cntrl file, set exptype
%%
% load cntrl
exptype='FF';
%%

fitmethod = 'T';
peranimal = 0;% if 1: calculates tau per animal as well as the combied one.
saveres = 0; % if 1: automatically makes timestamped folder and saves figures and tau estimates.

savdir = '/mnt/data/Mitra/figs/P2_L/bothdircombined/2020_02_24/LDA/version2/cosinesimilaritymats_and_timeconstants/saveres';
if strcmp(exptype,'FB')
    allanimals = [1:6];
    numanimals = 6;
elseif strcmp(exptype,'FF')
    allanimals = [1:7];
    numanimals = 7;
end


if strcmp(fitmethod,'T')
    fo = fitoptions('Method','NonlinearLeastSquares','Lower',[-inf,-inf,5],'Upper',[inf,inf,inf],'StartPoint',[0.5,1,1000],'MaxIter',1000);    
    % 5 was 0
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


% plotting params
scAlpha=0.5;
xjitter=0.1;
sc = 0; % whether to overlap scattered points

timescalems = 65; % used just for labeling xaxis -
nlags = 8;
binsizems = '150';
allanimalsn0cormat_go = [];
allanimalsn0cormat_nogo = [];
allanimalsn1cormat_go = [];
allanimalsn1cormat_nogo = [];
allanimalscomcormat_go = [];
allanimalscomcormat_nogo = [];
go_com_tc = nan(numanimals*nlags,nlags);
nogo_com_tc = nan(numanimals*nlags,nlags);
a_go_com_tc = nan(nlags,nlags);
a_nogo_com_tc = nan(nlags,nlags);
  
animaltau_com_go = {};
animaltau_com_nogo = {};
animaltau_com_c = {};

tau_com_go = {};
tau_com_nogo = {};
tau_com_c = {};

go_currentrow = 0;
nogo_currentrow = 0;


for animalnum = allanimals
    
    animal_go_currentrow = 0;
    animal_nogo_currentrow = 0;
    
    % go:
    Comcormat = nan(8,8);
    for i=1:numel(cntrl)
        Comcormat(:,:,end+1) = cntrl{i}.go{animalnum};
    end
    Comcormat(:,:,1)=[];
    Comcormat = nanmean(Comcormat,3);
    allanimalscomcormat_go = cat(3,allanimalscomcormat_go,Comcormat);    
    for ii=1:nlags
        go_currentrow = go_currentrow + 1;
        animal_go_currentrow = animal_go_currentrow + 1;
        for jj =1:nlags
            if (jj-ii) >= 0
                go_com_tc(go_currentrow,jj-ii+1) = Comcormat(ii,jj);
                a_go_com_tc(animal_go_currentrow,jj-ii+1) = Comcormat(ii,jj);
                
            end
        end
    end
 
    % nogo:   
    Comcormat = nan(8,8);
    for i=1:numel(cntrl)
        Comcormat(:,:,end+1) = cntrl{i}.nogo{animalnum};
    end
    Comcormat(:,:,1)=[];
    Comcormat = nanmean(Comcormat,3);
    allanimalscomcormat_nogo = cat(3,allanimalscomcormat_nogo,Comcormat);
    
    for ii=1:nlags
        nogo_currentrow = nogo_currentrow + 1;
        animal_nogo_currentrow = animal_nogo_currentrow + 1;
        for jj =1:nlags
            if (jj-ii) >= 0
                nogo_com_tc(nogo_currentrow,jj-ii+1) = Comcormat(ii,jj);
                a_nogo_com_tc(animal_nogo_currentrow,jj-ii+1) = Comcormat(ii,jj);
            end
        end
    end
    
    if peranimal
        figure;[animaltau_com_go{end+1},~]=X_cor_timeconstant_plot_helper(a_go_com_tc,'g',scAlpha,xjitter,sc,ftype,fo,timescalems); %untrained: (go_com_tc+nogo_com_tc)/2,'k'
        hold on;[animaltau_com_nogo{end+1},~]=X_cor_timeconstant_plot_helper(a_nogo_com_tc,'r',scAlpha,-xjitter,sc,ftype,fo,timescalems);
        title(['animal',num2str(animalnum),exptype,'-',binsizems,'ms'])
        
        figure;
        hold on;[animaltau_com_c{end+1},~]=X_cor_timeconstant_plot_helper((a_go_com_tc+a_nogo_com_tc)/2,'k',scAlpha,xjitter,sc,ftype,fo,timescalems); %untrained: (go_com_tc+nogo_com_tc)/2,'k'
        title(['animal-','combined-',num2str(animalnum),exptype,'-',binsizems,'ms'])
    end
    
    
end

%%% plotts:
% 2 styles of plots in helper: 1- data points, mean, and mean ci - you see the difference
%                              2- exponential fits with 0 intercept. shaded
%                              bars are prediction intervals of the
%                              function, fitting all simultaneously


figure;hold on;[tau6,~]=X_cor_timeconstant_plot_helper(nogo_com_tc,'r',scAlpha,-xjitter,sc,ftype,fo,timescalems,'exp');
hold on;[tau5,~]=X_cor_timeconstant_plot_helper(go_com_tc,'g',scAlpha,xjitter,sc,ftype,fo,timescalems,'exp'); 
title([exptype,'-',binsizems,'ms'])
tau_com_go{end+1} = tau5;
tau_com_nogo{end+1} = tau6;

figure;
hold on;[tau7,~]=X_cor_timeconstant_plot_helper((go_com_tc+nogo_com_tc)/2,'k',scAlpha,xjitter,sc,ftype,fo,timescalems,'exp'); 
title(['combined-',exptype,'-',binsizems,'ms'])
tau_com_c{end+1} = tau7;

%%%%% com cormats

figure;
s1=subplot(1,2,1);a1=imagesc(nanmean(allanimalscomcormat_go,3));
s2=subplot(1,2,2);a2=imagesc(nanmean(allanimalscomcormat_nogo,3));
s1.Title.String = (sprintf('go,%s,N = %d animals',exptype,size(allanimalscomcormat_go,3)));
s1.Colormap = colormap(jet(64));% redblue
set(s1,'CLim',[-0.7 0.7])
colorbar(a1.Parent)
s2.Title.String = (sprintf('nogo,%s,N = %d animals',exptype,size(allanimalscomcormat_go,3)));
s2.Colormap = colormap(jet(64));% redblue
set(s2,'CLim',[-0.7 0.7])
colorbar(a2.Parent)



if saveres
    cd(savdir)
    foldername = ['simul_',exptype,binsizems,'_',datestr(now,'yyyy_mm_dd_HH_MM_SS')];

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
        'animaltau_com_go','animaltau_com_nogo','animaltau_com_c',...
        'tau_com_go','tau_com_nogo','tau_com_c')
end