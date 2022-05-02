% This script combines V1 ceslls from ff and fb experiments. and shows
% average zscore of activitiy for different activity levels. This is to
% find a suitable threshold for setting activity to nan.

ff = load(fullfile('/mnt/data/Mitra/figs/P2_L/bothdircombined',['cells_','FF','_pl20_an150_lw20_exG1_onlyC1_onlyS0_plstyle1.mat']));
fb = load(fullfile('/mnt/data/Mitra/figs/P2_L/bothdircombined',['cells_','FB','_pl20_an150_lw20_exG1_onlyC1_onlyS0_plstyle1.mat']));
V1cells = [ff.V1cells,fb.V1cells];
LMcells = [ff.LMcells,fb.LMcells];
cell_keep_ind = nan;
[cell_keep_ind,LMcells,V1cells,Behcells,params] = removecells(cell_keep_ind,cellremovecriteria,LMcells,V1cells,Behcells,params);

bs = [];
bsstd = [];
gbs = cell2mat(cellfun(@(y) cellfun(@(x) nanmean(x),y.laAbs.go),V1cells,'UniformOutput',0));
ngbs = cell2mat(cellfun(@(y) cellfun(@(x) nanmean(x),y.laAbs.nogo),V1cells,'UniformOutput',0));
bs = [gbs/0.15,ngbs/0.15];
gbsstd = cell2mat(cellfun(@(y) cellfun(@(x) nanstd(x)/sqrt(length(x)),y.laAbs.go),V1cells,'UniformOutput',0));
ngbsstd = cell2mat(cellfun(@(y) cellfun(@(x) nanstd(x)/sqrt(length(x)),y.laAbs.nogo),V1cells,'UniformOutput',0));
bsstd = [gbsstd/0.15,ngbsstd/0.15];

% 
figure;hist(bs,1000)
figure;cdfplot(bsstd)

prctile(bs,20)
%% only V1, no LM cells
figure;scatter(bs,bs-2*bsstd,'.')
nbins = 100;
binsize = max(bs)/nbins;
val = nan(1,nbins);
for i=1:nbins
    ind=find(bs>((i-1)*binsize) & bs<((i)*binsize));
   % val(i) = nanmean(bs(ind)-2*bsstd(ind));
    val(i) = nanmean(bs(ind)/bsstd(ind));
end

figure;stairs(binsize*(1:nbins),val)
xlabel('baseline firing rate in the 150 ms bin');
ylabel('zscored activity averaged over cells and bins');
hold on;line([0 100],[2,2],'Color','k')