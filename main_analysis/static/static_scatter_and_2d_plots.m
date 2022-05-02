% load cells, run the first 2 sections of staticpropertiesSilencing.m after setting target ( = remove nans and inf as well)
% then set th following 2 params
% scatters with groups of depth or retin , a pre-set number of cells are
% selected (to make groups comparable)
%% params
exptype = 'FF';
retintype = 'shank';% 'monitor', 'shank'(value of allretinlocald is either distance to local shanks or 
% long range (center of silencing) depending on staticpropertiesVssilencing.m)

if strcmp(retintype,'shank')
    allcenterdistance = allretinlocald; 
end
%%
figure
count = 0;
cmap = jet(4);
for i=min(alldepth):200:max(alldepth)

    ind=find(alldepth>=i & alldepth<i+100);
    if numel(ind)
            count = count+1;
    ind = ind(1:min(30,numel(ind)));
     hold on;scatter([allgobs(ind) allnogobs(ind)],[allgo(ind) allnogo(ind)],100,'.','MarkerEdgeColor',cmap(count,:),...
     'MarkerFaceColor','none')
    end
end
title('depth')


figure

count = 0;
cmap = jet(8);
for i=0:.5:3

    ind=find(allcenterdistance>=i & allcenterdistance<i+1);
    if numel(ind)
            count = count+1;
            ind = ind(1:min(10,numel(ind)));
     hold on;scatter([allgobs(ind) allnogobs(ind)],[allgo(ind) allnogo(ind)],100,'.','MarkerEdgeColor',cmap(count,:),...
     'MarkerFaceColor','none')
    end
end
title('rf-center-distance')

%% 2d plots : load cells, run the first 2 sections of staticpropertiesSilencing.m (remove nans and inf as wellß)
figure;
nx = 10;
ny = 5;
%%% assuming 150 ms an bin
binsizes = 150/1000; % all firing rates get converted to Hz
plotd = 2; % 2:2d, 1:1d(average over baseline)
%%%%%%%%%%%%%% depth

xvec = [allgobs allnogobs]/binsizes;
yvec = [alldepth,alldepth];
zvec_ls = [allgo allnogo]/binsizes;
zvec_delta = ([allgo allnogo] - [allgobs allnogobs])/binsizes; 
zvec_pc = 100*([allgo allnogo] - [allgobs allnogobs])./[allgobs allnogobs];


xmesh=linspace(0,7,nx)/binsizes;%
ymesh =  linspace(min(yvec),max(yvec),ny);
xstep=unique(diff(xmesh));
ystep=unique(diff(ymesh));
if numel(xstep)>1
    xstep = xstep(1);
end
if numel(ystep)>1
    ystep = ystep(1);
end

p_ls = nan(nx,ny);
p_delta = nan(nx,ny);
p_pc = nan(nx,ny);

for i = 1:length(xmesh)
    for j = 1:length(ymesh)
        xind = find(xvec>=(xmesh(i)-xstep/2)&(xvec<(xmesh(i)+xstep/2)));
        yind = find(yvec>=(ymesh(j)-ystep/2)&(yvec<(ymesh(j)+ystep/2)));
        p_ls(i,j) = nanmean(zvec_ls(intersect(intersect(xind,yind),find(~isinf(zvec_ls)))));
        p_delta(i,j) = nanmean(zvec_delta(intersect(intersect(xind,yind),find(~isinf(zvec_delta)))));
        p_pc(i,j) = nanmean(zvec_pc(intersect(intersect(xind,yind),find(~isinf(zvec_pc)))));   
    end
end



s=subplot(2,3,1);
if plotd ==2 
    imagesc(xmesh,ymesh,p_ls')
    colorbar
    ax = gca;
    if strcmp(exptype,'FF')
        ax.CLim = [0 32];
    elseif strcmp(exptype,'FB')
        ax.CLim = [0 41];
    end
else
    plot(nanmean(p_ls,1),ymesh,'k');set(gca,'Ydir','reverse')
end
s.Title.String = ' silencing activity';
xlabel('basline activity')
ylabel('depth')


s=subplot(2,3,2);
if plotd == 2
    imagesc(xmesh,ymesh,p_delta')
    colorbar
else
    plot(nanmean(p_delta,1),ymesh,'k');set(gca,'Ydir','reverse')
end

s.Title.String = ' delta activity';
xlabel('basline activity')
ylabel('depth')


s=subplot(2,3,3);
if plotd == 2
    imagesc(xmesh,ymesh,p_pc')
    colorbar
else
    plot(nanmean(p_pc,1),ymesh,'k');set(gca,'Ydir','reverse')
end
s.Title.String = ' percentage change';
xlabel('basline activity')
ylabel('depth')

colormap([1 1 1;parula(64)]);
%%%%%%%%%%%%%% retin
p_ls = nan(nx,ny);
p_delta = nan(nx,ny);
p_pc = nan(nx,ny);

yvec = [allcenterdistance,allcenterdistance];
ymesh = linspace(min(yvec),2.2,ny);
ystep=unique(diff(ymesh));
if numel(ystep)>1
    ystep = ystep(1);
end

for i = 1:length(xmesh)
    for j = 1:length(ymesh)
        xind = find(xvec>=(xmesh(i)-xstep/2)&(xvec<(xmesh(i)+xstep/2)));
        yind = find(yvec>=(ymesh(j)-ystep/2)&(yvec<(ymesh(j)+ystep/2)));
        p_ls(i,j) = nanmean(zvec_ls(intersect(intersect(xind,yind),find(~isinf(zvec_ls)))));
        p_delta(i,j) = nanmean(zvec_delta(intersect(intersect(xind,yind),find(~isinf(zvec_delta)))));        
        p_pc(i,j) = nanmean(zvec_pc(intersect(intersect(xind,yind),find(~isinf(zvec_pc)))));
    end
end



s=subplot(2,3,4);
if plotd ==2 
    imagesc(xmesh,ymesh,p_ls')
    colorbar
    ax = gca;
    if strcmp(exptype,'FF')
        ax.CLim = [0 32];
    elseif strcmp(exptype,'FB')
        ax.CLim = [0 41];
    end
else
    plot(nanmean(p_ls,1),ymesh,'k');set(gca,'Ydir','reverse')
end
s.Title.String = ' silencing activity';
xlabel('basline activity')
ylabel('receptive field center distance')

s=subplot(2,3,5);
if plotd == 2
    imagesc(xmesh,ymesh,p_delta')
    colorbar
else
    plot(nanmean(p_delta,1),ymesh,'k');set(gca,'Ydir','reverse')
end
s.Title.String = ' delta activity';
xlabel('basline activity')
ylabel('receptive field center distance')

s=subplot(2,3,6);
if plotd == 2
    imagesc(xmesh,ymesh,p_pc')
    colorbar
else
    plot(nanmean(p_pc,1),ymesh,'k');set(gca,'Ydir','reverse')
end
s.Title.String = ' percentage change';
xlabel('basline activity')
ylabel('receptive field center distance')

%colormap([parula(64)]);
colormap([1 1 1;parula(64)]);

%% relation between depth and retin themselves: for one condition 

figure;scatter(allcenterdistance,alldepth,100,'k.')
[r_rd,p_rd] = corrcoef(allcenterdistance,alldepth);
title(sprintf('r-sq=%f, p=%d',r_rd(1,2)^2,p_rd(1,2)))
ylabel('depth,- up');xlabel('rf center distance from center - 1~20 deg')
%% show bs vs ls relation
binsize = 150/1000;
figure;scatter([allgobs,allnogobs]/binsize,[allgo,allnogo]/binsize,50,'k.');
xlim([0 70]);ylim([0 70]); xlabel('bs-Hz');ylabel('ls-Hz')

figure;scatter([alldepth alldepth],([allgo,allnogo]-[allgobs,allnogobs])./([allgobs,allnogobs]),50,'k.')

%% extra plots
binsize = 150/1000;
figure;scatter([allgobs,allnogobs]/binsize,100*([allgo,allnogo]-[allgobs,allnogobs])./([allgobs,allnogobs]),'k.');
xlabel('baseline firing rate (hz)')
ylabel('silencing effect(%)')

figure;
scatter([alldepth alldepth],100*([allgo,allnogo]-[allgobs,allnogobs])./([allgobs,allnogobs]),50,'k.')
effect = 100*([allgo,allnogo]-[allgobs,allnogobs])./([allgobs,allnogobs]);
%effect = [allgobs,allnogobs]/binsize;
effect(find(isinf(effect)))= nan;
dep = [alldepth alldepth];
figure;
valeffect = [];
valbas = [];
for d = unique(alldepth)
    in = find(dep == d);
    hold on;errorbar(d,nanmedian(effect(in)),nanmedian(effect(in))-prctile(effect(in),2.5),...
        prctile(effect(in),97.5)-nanmedian(effect(in)),'k')
    hold on;scatter(d,nanmean(effect(in)),'.k')
    valeffect(end+1)=nanmean(effect(in));
    
end
hold on;plot(unique(alldepth),valeffect,'k')
xlabel('depth')
ylabel('silencing effect(%)')

