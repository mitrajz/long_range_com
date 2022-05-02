% load cells, run the first 2-3 sections of staticBehaviorpropertiesSilencing.m after setting target ( = remove nans and inf as well)
% then set th followingparams
% scatters with groups of depth or retin , a pre-set number of cells are
% selected (to make groups comparable)
%% params
exptype = 'FB';
%% 2d plots : load cells, run the first 2 sections of staticpropertiesSilencing.m (remove nans and inf as wellß)
figure;
nx = 10;
ny =6;
%%% assuming 150 ms an bin
binsizes = 150/1000; % all firing rates get converted to Hz
plotd = 2; % 2:2d, 1:1d(average over baseline)
%%%%%%%%%%%%%% depth

xvec = [allgobs allnogobs]/binsizes;
yvec = [all_earlystim_tun,all_earlystim_tun];
zvec_ls = [allgo allnogo]/binsizes;
zvec_delta = ([allgo allnogo] - [allgobs allnogobs])/binsizes; 
zvec_pc = 100*([allgo allnogo] - [allgobs allnogobs])./[allgobs allnogobs];


xmesh=linspace(0,7,nx)/binsizes;
ymesh =  linspace(0.1,0.7,ny);% change to 0.1 if you have to do earlystim
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



s=subplot(2,1,1);
if plotd ==2 
    imagesc(xmesh,ymesh,p_ls')
    colorbar
    ax = gca;
    if strcmp(exptype,'FF')
        ax.CLim = [0 32];
    elseif strcmp(exptype,'FB')
        ax.CLim = [0 48];
    end
else
    plot(nanmean(p_ls,1),ymesh,'k');set(gca,'Ydir','reverse')
end
s.Title.String = ' silencing activity';
xlabel('basline activity')
ylabel('all_earlystim_tun')
colormap([1 1 1;parula(64)]);
%%%% retin


p_ls = nan(nx,ny);
p_delta = nan(nx,ny);
p_pc = nan(nx,ny);

yvec = [all_lick_tun,all_lick_tun];
%ymesh = linspace(min(yvec),max(yvec),ny);%linspace(min(yvec),max(yvec),ny);
ymesh =  linspace(-0.5,0.5,ny);
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



s=subplot(2,1,2);
if plotd ==2 
    imagesc(xmesh,ymesh,p_ls')
    colorbar
    ax = gca;
    if strcmp(exptype,'FF')
        ax.CLim = [0 28];
    elseif strcmp(exptype,'FB')
        ax.CLim = [0 48];
    end
else
    plot(nanmean(p_ls,1),ymesh,'k');set(gca,'Ydir','reverse')
end
s.Title.String = ' silencing activity';
xlabel('basline activity')
ylabel('all_lick_tun')

colormap([1 1 1;parula(64)]);

