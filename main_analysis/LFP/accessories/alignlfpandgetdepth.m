function   [shankdepthmap,refdepth,lfpdepth_f] = alignlfpandgetdepth(order,data_lfp_bs,animalshankorder,lfptimewindow,smoothinglength);
% note that lfp here is broad band (<300Hz), filtering it doesn't make a
% considerable difference in csd
refdepth = cell(1,numel(order));
lfpdepth_f = cell(1,numel(order));
shankdepthmap = cell(1,numel(order));
for shanknum = 1:numel(order)
    % just the go trials:
    LAvchdata = squeeze(nanmean(data_lfp_bs.go(find(strcmp(animalshankorder,order{shanknum})),:,:),2));
    figure;imagesc((lfptimewindow(1):lfptimewindow(2))/30,[],LAvchdata);
    dedchs = input('which channels are dead? ([x y], if none [])');
    LAvchdata = fillindedchanels(LAvchdata,dedchs);
    clear SLAvchdata1 SLAvchdata2
    % smooth or median or triangular kernel? here 2 methods are used
    for i=1:size(LAvchdata,2)
        SLAvchdata1(:,i)=smooth(LAvchdata(:,i),smoothinglength);
        SLAvchdata2(:,i)=medfilt1(LAvchdata(:,i),smoothinglength,'truncate');
    end
    lfpdepth_f{shanknum}=figure;
    lfpdepth_f{shanknum}.Units= 'normalized';
    lfpdepth_f{shanknum}.Position = [ 0 0.5 1 0.5];
    lfpdepth_f{shanknum}.Name = order{shanknum};
    subplot(1,3,1);si1=imagesc((lfptimewindow(1):lfptimewindow(2))/30,[],SLAvchdata1);
    xlim([-50 200])
    subplot(1,3,2);si2=imagesc((lfptimewindow(1):lfptimewindow(2))/30,[],-diff(SLAvchdata1,2,1));
    xlim([-50 200])
    si2.Parent.CLim = [-50 50];
    subplot(1,3,3);si3=imagesc((lfptimewindow(1):lfptimewindow(2))/30,[],-diff(SLAvchdata2,2,1));
    xlim([-50 200])
    si3.Parent.CLim = [-50 50];
    %%%%%%%%%%%
    si1.Parent.Title.String = sprintf('%s - smoothlfp - go',order{shanknum});
    si1.Parent.XLabel.String = 'time (ms)';
    si1.Parent.YLabel.String = 'ch number';
    si2.Parent.Title.String = sprintf('%s - csd(smooth) - go',order{shanknum});
    si2.Parent.XLabel.String = 'time (ms)';
    si2.Parent.YLabel.String = 'depth from L4 base';
    si3.Parent.Title.String = sprintf('%s - csd(med) - go',order{shanknum});
    si3.Parent.XLabel.String = 'time (ms)';
    si3.Parent.YLabel.String = 'depth from L4 base';
    %%%%%%%%%%%
    refdepth{shanknum} = input('where is the base of L4?');
    shankdepthmap{shanknum} = 25*((1:numel(find(strcmp(animalshankorder,order{shanknum})))) -refdepth{shanknum}); % 25 micros
    si2.Parent.YTick=1:1:numel(find(strcmp(animalshankorder,order{shanknum})))-2;
    si2.Parent.YTickLabel=shankdepthmap{shanknum}(1:end-2);
    hold(si2.Parent,'on');line(si2.Parent,[-50 200] , [refdepth{shanknum}, refdepth{shanknum}],'Color','k')
    hold(si2.Parent,'on');line(si2.Parent,[-50 200] , [refdepth{shanknum}-4, refdepth{shanknum}-4],'Color','k') % 100 microns above: roughly top of L4
    si3.Parent.YTick=1:1:numel(find(strcmp(animalshankorder,order{shanknum})))-2;
    si3.Parent.YTickLabel=shankdepthmap{shanknum}(1:end-2);
    hold(si3.Parent,'on');line(si3.Parent,[-50 200] , [refdepth{shanknum}, refdepth{shanknum}],'Color','k')
    hold(si3.Parent,'on');line(si3.Parent,[-50 200] , [refdepth{shanknum}-4, refdepth{shanknum}-4],'Color','k') % 100 microns above: roughly top of L$
    
end
end