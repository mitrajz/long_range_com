%% This code is currently not used in any figures/sups

mean_standarderr_coeff =2.73%1.96; 
nrep = 100;% 100
%% first load cells , only clean (according to latest methods) R0C1
exptype = 'FB';
load(['cells_',exptype,'_pl20_an150_lw20_exG1_onlyC1_onlyS0_plstyle1.mat']);


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
cellremovecriteria.sigTH = 0.05;
cellremovecriteria.minspkTH = 1; %responsiveness param (method 4 only)- def:1
cellremovecriteria.showexampleplots = 0; % order:targetgo, targetnogo, not_target go, not_targetnogo

[LMcells,V1cells,Behcells,params] = cleanupcells(cleanupcriteria,LMcells,V1cells,Behcells,params);
[cell_keep_ind,LMcells,V1cells,Behcells,params] = removecells(cell_keep_ind,cellremovecriteria,LMcells,V1cells,Behcells,params);

%% set target cells and assess significance (per cell, thereforeBonferoni correction)
if strcmp(exptype,'FF')
    targetcell = LMcells;
elseif strcmp(exptype,'FB')
    targetcell = V1cells;
end

allgo = cell2mat(cellfun(@(x) x.smb.go, targetcell,'UniformOutput',0));
allnogo = cell2mat(cellfun(@(x) x.smb.nogo, targetcell,'UniformOutput',0));
all = cell2mat([cellfun(@(x) x.smb.nogo, targetcell,'UniformOutput',0),cellfun(@(x) x.smb.go, targetcell,'UniformOutput',0)]);

% standard error of the mean
standarderr_prct_go = nan(length(targetcell),length(targetcell{1}.laAbs.go));
standarderr_prct_nogo = nan(length(targetcell),length(targetcell{1}.laAbs.nogo));
for i=1:length(targetcell)
    % go
    for l=1:length(targetcell{i}.laAbs.go)
        mean_standarderr = nanstd(targetcell{i}.laAbs.go{l})/sqrt(numel(targetcell{i}.laAbs.go{l}));
        standarderr_prct_go(i,l) = 100 * mean_standarderr_coeff * mean_standarderr/nanmean(targetcell{i}.laAbs.go{l});
        if isinf(standarderr_prct_go(i,l))
            standarderr_prct_go(i,l) = nan;
        end
    end
    % nogo
    for l=1:length(targetcell{i}.laAbs.nogo)
        mean_standarderr = nanstd(targetcell{i}.laAbs.nogo{l})/sqrt(numel(targetcell{i}.laAbs.nogo{l}));
        standarderr_prct_nogo(i,l) = 100 * mean_standarderr_coeff * mean_standarderr/nanmean(targetcell{i}.laAbs.nogo{l});
        if isinf(standarderr_prct_nogo(i,l))
            standarderr_prct_nogo(i,l) = nan;
        end
    end
end

allmeanconfs = [reshape(standarderr_prct_nogo',1,[]),reshape(standarderr_prct_go',1,[])];
gomeanconfs = reshape(standarderr_prct_go',1,[]);
nogomeanconfs = reshape(standarderr_prct_nogo',1,[]);

all_s = all; allgo_s = allgo; allnogo_s = allnogo;
all_s(abs(all_s)<allmeanconfs) = nan;
allgo_s(abs(allgo_s)<gomeanconfs) = nan;
allnogo_s(abs(allnogo_s)<nogomeanconfs) = nan;



%% pie charts
% this is how many cells were affected at least at one point.
% color correct the pie charts across ff and fb posthoc

nsigbins_go = [numel(find(sum(isnan(reshape(allgo_s,8,[])),1)==8))/numel(targetcell),...
    numel(find(sum(isnan(reshape(allgo_s,8,[])),1)==7))/numel(targetcell),...
    numel(find(sum(isnan(reshape(allgo_s,8,[])),1)==6))/numel(targetcell),...
    numel(find(sum(isnan(reshape(allgo_s,8,[])),1)==5))/numel(targetcell),...
    numel(find(sum(isnan(reshape(allgo_s,8,[])),1)==4))/numel(targetcell),...
    numel(find(sum(isnan(reshape(allgo_s,8,[])),1)==3))/numel(targetcell),...
    numel(find(sum(isnan(reshape(allgo_s,8,[])),1)==2))/numel(targetcell),...
    numel(find(sum(isnan(reshape(allgo_s,8,[])),1)==1))/numel(targetcell),...
    numel(find(sum(isnan(reshape(allgo_s,8,[])),1)==0))/numel(targetcell)];
atleast1_go = numel(find(sum(isnan(reshape(allgo_s,8,[])),1)<8))/numel(targetcell);

nsigbins_nogo = [numel(find(sum(isnan(reshape(allnogo_s,8,[])),1)==8))/numel(targetcell),...
    numel(find(sum(isnan(reshape(allnogo_s,8,[])),1)==7))/numel(targetcell),...
    numel(find(sum(isnan(reshape(allnogo_s,8,[])),1)==6))/numel(targetcell),...
    numel(find(sum(isnan(reshape(allnogo_s,8,[])),1)==5))/numel(targetcell),...
    numel(find(sum(isnan(reshape(allnogo_s,8,[])),1)==4))/numel(targetcell),...
    numel(find(sum(isnan(reshape(allnogo_s,8,[])),1)==3))/numel(targetcell),...
    numel(find(sum(isnan(reshape(allnogo_s,8,[])),1)==2))/numel(targetcell),...
    numel(find(sum(isnan(reshape(allnogo_s,8,[])),1)==1))/numel(targetcell),...
    numel(find(sum(isnan(reshape(allnogo_s,8,[])),1)==0))/numel(targetcell)];
atleast1_nogo = numel(find(sum(isnan(reshape(allnogo_s,8,[])),1)<8))/numel(targetcell);

nsigbins_av = nanmean([nsigbins_go;nsigbins_nogo]);

f=figure;subplot(1,3,1);pie(nsigbins_go);
ax1 = gca;ax1.Colormap = [0.5,0.5,0.5;parula(8)];title(ax1,['go-',num2str(atleast1_go),'at leat 1 bin sig']);
legend(ax1,{'0','1','2','3','4','5','6','7','8'})
subplot(1,3,2);pie(nsigbins_nogo);
ax2 = gca;ax2.Colormap = [0.5,0.5,0.5;parula(8)];title(ax2,['nogo-',num2str(atleast1_nogo),'at leat 1 bin sig']);
legend(ax2,{'0','1','2','3','4','5','6','7','8'})
subplot(1,3,3);pie(nsigbins_av);
ax3 = gca;ax3.Colormap = [0.5,0.5,0.5;parula(8)];title(ax3,['averaged-',num2str(nanmean([atleast1_go,atleast1_nogo])),'at leat 1 bin sig']);
legend(ax3,{'0','1','2','3','4','5','6','7','8'})

%% see which bins were significant
figure;plot(nansum([sum(~isnan(reshape(allgo_s,8,[])),2),sum(~isnan(reshape(allnogo_s,8,[])),2)],2));ylim([0 200])

%% signal correlation between 2 halves of influence traces over time, compared to random chance

histbins = -1:0.15:1;% 0.15
r_histbins = -1:0.15:1; % smoother bins?
cellsel = '1sig';%'all' or '1sig'
shufflenrep = 10;%100
shuffletype = 'all'; % 'pertrace' or 'all'
% go
if strcmp(cellsel,'1sig')
    cellindcs = find(sum(isnan(reshape(allgo_s,8,[])),1)<8);
else
    cellindcs = 1:numel(allgo_s)/8;
end
bs= cellfun(@(x) x.laAbs.go,targetcell(cellindcs),'UniformOutput',0);
ls= cellfun(@(x) x.laAls.go,targetcell(cellindcs),'UniformOutput',0);
[self_sim_go] = calculateselfsim(bs,ls,0,shufflenrep,shuffletype);
[r_self_sim_go] = calculateselfsim(bs,ls,1,shufflenrep,shuffletype);

% nogo
if strcmp(cellsel,'1sig')
    cellindcs = find(sum(isnan(reshape(allnogo_s,8,[])),1)<8);
else
    cellindcs = 1:numel(allnogo_s)/8;
end
bs= cellfun(@(x) x.laAbs.nogo,targetcell(cellindcs),'UniformOutput',0);
ls= cellfun(@(x) x.laAls.nogo,targetcell(cellindcs),'UniformOutput',0);
[self_sim_nogo] = calculateselfsim(bs,ls,0,shufflenrep,shuffletype);
[r_self_sim_nogo] = calculateselfsim(bs,ls,1,shufflenrep,shuffletype);


self_sim_pooled = [self_sim_go,self_sim_nogo];
r_self_sim_pooled = [r_self_sim_go,r_self_sim_nogo];

%%% plots
figure;subplot(1,3,1);
hold on;histogram(r_self_sim_go,r_histbins,'Normalization','pdf','FaceColor','k','EdgeAlpha',0);
hold on;histogram(self_sim_go,histbins,'Normalization','pdf','FaceColor','g','EdgeAlpha',0)
hold on;line([nanmedian(self_sim_go),nanmedian(self_sim_go)],[0,3],'Color','g')
hold on;line([nanmedian(r_self_sim_go),nanmedian(r_self_sim_go)],[0,3],'Color','k')
hold on;text(nanmedian(self_sim_go),3,num2str(nanmedian(self_sim_go)))
hold on;text(nanmedian(r_self_sim_go),2.8,num2str(nanmedian(r_self_sim_go)))
hold on;text(-1,2.5,['ranksum p=', num2str(ranksum(self_sim_go,r_self_sim_go))])

subplot(1,3,2);
hold on;histogram(r_self_sim_nogo,r_histbins,'Normalization','pdf','FaceColor','k','EdgeAlpha',0);
hold on;histogram(self_sim_nogo,histbins,'Normalization','pdf','FaceColor','r','EdgeAlpha',0);
hold on;line([nanmedian(self_sim_nogo),nanmedian(self_sim_nogo)],[0,3],'Color','r')
hold on;line([nanmedian(r_self_sim_nogo),nanmedian(r_self_sim_nogo)],[0,3],'Color','k')
hold on;text(nanmedian(self_sim_nogo),3,num2str(nanmedian(self_sim_nogo)))
hold on;text(nanmedian(r_self_sim_nogo),2.8,num2str(nanmedian(r_self_sim_nogo)))
hold on;text(-1,2.5,['ranksum p=', num2str(ranksum(self_sim_nogo,r_self_sim_nogo))])



subplot(1,3,3);
hold on;histogram(r_self_sim_pooled,r_histbins,'Normalization','pdf','FaceColor','k','EdgeAlpha',0);
hold on;histogram(self_sim_pooled,histbins,'Normalization','pdf','FaceColor','b','EdgeAlpha',0);
hold on;line([nanmedian(self_sim_pooled),nanmedian(self_sim_pooled)],[0,3],'Color','b')
hold on;line([nanmedian(r_self_sim_pooled),nanmedian(r_self_sim_pooled)],[0,3],'Color','k')
hold on;text(nanmedian(self_sim_pooled),3,num2str(nanmedian(self_sim_pooled)))
hold on;text(nanmedian(r_self_sim_pooled),2.8,num2str(nanmedian(r_self_sim_pooled)))
hold on;text(-1,2.5,['ranksum p=', num2str(ranksum(self_sim_pooled,r_self_sim_pooled))])

%% cross sig cors

% go
if strcmp(cellsel,'1sig')
    cellindcs = find(sum(isnan(reshape(allgo_s,8,[])),1)<8);
else
    cellindcs = 1:numel(allgo_s)/8;
end
bs= cellfun(@(x) x.laAbs.go,targetcell(cellindcs),'UniformOutput',0);
ls= cellfun(@(x) x.laAls.go,targetcell(cellindcs),'UniformOutput',0);
[x_sim_go] = calculatecrossim(bs,ls,0,shufflenrep,shuffletype);
[r_x_sim_go] = calculatecrossim(bs,ls,1,shufflenrep,shuffletype);

% nogo
if strcmp(cellsel,'1sig')
    cellindcs = find(sum(isnan(reshape(allnogo_s,8,[])),1)<8);
else
    cellindcs = 1:numel(allnogo_s)/8;
end
bs= cellfun(@(x) x.laAbs.nogo,targetcell(cellindcs),'UniformOutput',0);
ls= cellfun(@(x) x.laAls.nogo,targetcell(cellindcs),'UniformOutput',0);
[x_sim_nogo] = calculatecrossim(bs,ls,0,shufflenrep,shuffletype);
[r_x_sim_nogo] = calculatecrossim(bs,ls,1,shufflenrep,shuffletype);


x_sim_pooled = [x_sim_go,x_sim_nogo];
r_x_sim_pooled = [r_x_sim_go,r_x_sim_nogo];

%%% plots
figure;subplot(1,3,1);
hold on;histogram(r_x_sim_go,r_histbins,'Normalization','pdf','FaceColor','k','EdgeAlpha',0);
hold on;histogram(x_sim_go,histbins,'Normalization','pdf','FaceColor','g','EdgeAlpha',0)
hold on;line([nanmedian(x_sim_go),nanmedian(x_sim_go)],[0,3],'Color','g')
hold on;line([nanmedian(r_x_sim_go),nanmedian(r_x_sim_go)],[0,3],'Color','k')
hold on;text(nanmedian(x_sim_go),3,num2str(nanmedian(x_sim_go)))
hold on;text(nanmedian(r_x_sim_go),2.8,num2str(nanmedian(r_x_sim_go)))
hold on;text(-1,2.5,['ranksum p=', num2str(ranksum(x_sim_go,r_x_sim_go))])

subplot(1,3,2);
hold on;histogram(r_x_sim_nogo,r_histbins,'Normalization','pdf','FaceColor','k','EdgeAlpha',0);
hold on;histogram(x_sim_nogo,histbins,'Normalization','pdf','FaceColor','r','EdgeAlpha',0);
hold on;line([nanmedian(x_sim_nogo),nanmedian(x_sim_nogo)],[0,3],'Color','r')
hold on;line([nanmedian(r_x_sim_nogo),nanmedian(r_x_sim_nogo)],[0,3],'Color','k')
hold on;text(nanmedian(x_sim_nogo),3,num2str(nanmedian(x_sim_nogo)))
hold on;text(nanmedian(r_x_sim_nogo),2.8,num2str(nanmedian(r_x_sim_nogo)))
hold on;text(-1,2.5,['ranksum p=', num2str(ranksum(x_sim_nogo,r_x_sim_nogo))])



subplot(1,3,3);
hold on;histogram(r_x_sim_pooled,r_histbins,'Normalization','pdf','FaceColor','k','EdgeAlpha',0);
hold on;histogram(x_sim_pooled,histbins,'Normalization','pdf','FaceColor','b','EdgeAlpha',0);
hold on;line([nanmedian(x_sim_pooled),nanmedian(x_sim_pooled)],[0,3],'Color','b')
hold on;line([nanmedian(r_x_sim_pooled),nanmedian(r_x_sim_pooled)],[0,3],'Color','k')
hold on;text(nanmedian(x_sim_pooled),3,num2str(nanmedian(x_sim_pooled)))
hold on;text(nanmedian(r_x_sim_pooled),2.8,num2str(nanmedian(r_x_sim_pooled)))
hold on;text(-1,2.5,['ranksum p=', num2str(ranksum(x_sim_pooled,r_x_sim_pooled))])


%% only for 1 trial type, go or nogo
function [self_sim] = calculateselfsim(bs,ls,dornd,shufflenrep,shuffletype)



for celn = 1:numel(ls)
    mintr =  min(cellfun(@(x) numel(x),ls{celn}));
    for l = 1:8
        ls{celn}{l} = ...
            100*(ls{celn}{l}(1:mintr) - nanmean(bs{celn}{l}))/...
            nanmean(bs{celn}{l});
    end
end

if dornd
    nr = shufflenrep;
else
    nr = 1;
end

self_sim=nan(nr,numel(ls));



for perm = 1:nr
    for celn = 1:numel(ls)
        % ls activity mintr*8
        a = cell2mat(ls{celn});
        even = reshape(a(2:2:end,1:8)',1,[]);
        odd = reshape(a(1:2:end,1:8)',1,[]);
        
        if numel(even)>numel(odd)
            even(end-7:end) = [];
        elseif numel(odd)>numel(even)
            odd(end-7:end) = [];
        end
        
        % shuffling across the whole trace
        if dornd
            rng(perm,'twister');
            if strcmp(shuffletype,'all')
                even = even(randperm(numel(even)));
                odd = odd(randperm(numel(odd)));
            else
                even = reshape(a(2:2:end,randperm(8))',1,[]);
                odd = reshape(a(1:2:end,randperm(8))',1,[]);
                
                if numel(even)>numel(odd)
                    even(end-7:end) = [];
                elseif numel(odd)>numel(even)
                    odd(end-7:end) = [];
                end
            end
        end

        nanindcs = union(find(isnan(even)),find(isnan(odd)));
        odd(nanindcs) = [];
        even(nanindcs) = [];
        r = corrcoef(even,odd);
        try
            self_sim(perm,celn) = r(1,2);
        catch
            self_sim(perm,celn)= nan;
        end
    end
end

%self_sim = nanmean(self_sim,1);
self_sim = reshape(self_sim,1,[]);
end
%%
function [x_sim] = calculatecrossim(bs,ls,dornd,shufflenrep,shuffletype)

if dornd
    nr = shufflenrep;
else
    nr = 1;
end

x_sim=nan(nr,((numel(ls)^2-numel(ls))/2));


for celn = 1:numel(ls)
    mintr =  min(cellfun(@(x) numel(x),ls{celn}));
    for l = 1:8
        ls{celn}{l} = ...
            100*(ls{celn}{l}(1:mintr) - nanmean(bs{celn}{l}))/...
            nanmean(bs{celn}{l});
    end
end




for perm = 1:nr
    cellcount = 1;
    for celn = 1:numel(ls)
        for celm = celn+1:numel(ls)      
        % ls activity mintr*8
        a_n = cell2mat(ls{celn});
        even_n = reshape(a_n(2:2:end,1:8)',1,[]);
        odd_n = reshape(a_n(1:2:end,1:8)',1,[]);
        
        a_m = cell2mat(ls{celm});
        even_m = reshape(a_m(2:2:end,1:8)',1,[]);
        odd_m = reshape(a_m(1:2:end,1:8)',1,[]);
        
        trunc = min([numel(even_n),numel(odd_n),numel(even_m),numel(odd_m)]);
        even_m(trunc+1:end) = [];
        even_n(trunc+1:end) = [];
        odd_m(trunc+1:end) = [];
        odd_n(trunc+1:end) = [];
   
            % shuffling across the whole trace
            if dornd
                if strcmp(shuffletype,'all')
                    rng(perm,'twister');
                    even_n = even_n(randperm(numel(even_n)));
                    odd_n = odd_n(randperm(numel(odd_n)));
                    even_m = even_m(randperm(numel(even_m)));
                    odd_m = odd_m(randperm(numel(odd_m)));
                else
                    even_n = reshape(a_n(2:2:end,randperm(8))',1,[]);
                    odd_n = reshape(a_n(1:2:end,randperm(8))',1,[]);
                    even_m = reshape(a_m(2:2:end,randperm(8))',1,[]);
                    odd_m = reshape(a_m(1:2:end,randperm(8))',1,[]);
                    
                    trunc = min([numel(even_n),numel(odd_n),numel(even_m),numel(odd_m)]);
                    even_m(trunc+1:end) = [];
                    even_n(trunc+1:end) = [];
                    odd_m(trunc+1:end) = [];
                    odd_n(trunc+1:end) = [];
                    
                end
            end
        


        nanindcs = union(find(isnan(even_m)),find(isnan(odd_n)));
        odd_n(nanindcs) = [];
        even_m(nanindcs) = [];
        
        nanindcs = union(find(isnan(even_n)),find(isnan(odd_m)));
        odd_m(nanindcs) = [];
        even_n(nanindcs) = [];
        
        r1 = corrcoef(even_m,odd_n);
        r2 = corrcoef(even_n,odd_m);
        try
            x_sim(perm,cellcount) = (r1(1,2)+r2(1,2))/2;

        catch
            x_sim(perm,cellcount)= nan;
        end
        cellcount = cellcount + 1;
        end
    end
end

%x_sim = nanmean(x_sim,1);
x_sim = reshape(x_sim,1,[]);
end
