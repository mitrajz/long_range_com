
figure
all =[];
for i = 1:length(LMcells)
targetcell = LMcells{i};
    difs = cellfun(@(x,y) 100*(x- nanmean(y))./nanmean(y),targetcell.laAls.go,targetcell.laAbs.go,'UniformOutput',0);
all = [all;cellfun(@(x) nanmean(x),difs)];

end
all(isinf(all))=nan;
all(find(isnan(mean(all,2))),:)=[];
[~,ind] = sort(nanmean(all,2));
subplot(1,2,1);imagesc(all(ind,:));

[cmap]=cbrewer('div', 'RdBu',50);
colormap(flipud(cmap))
set(gca,'clim',[-100,100])
colorbar
title('ff-go')
subplot(1,2,2);imagesc((std(all(ind,:)')./1  )')
 set(gca,'Colormap',parula(50))
colorbar


figure; plot(nanmean(all(ind,:),2),(nanstd(all(ind,:)')))
xlabel('average effect');ylabel('standard deviation of effect')