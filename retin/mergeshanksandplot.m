function mergeshanksandplot(unit_maps_sh,unittype,h1,h2,areaname,Normalize)
% unittype: 'all', 'KSgood', 'taskmatched'
% Normalize: normalize each units map to its max or not


%%% merge 2 shanks: useful mainly for plotting, not saved

unit_maps=unit_maps_sh{1};
fields = fieldnames(unit_maps_sh{1});
for k=1:numel(fields)
    afield = fields{k};
    for n = 2:length(unit_maps_sh)
        unit_maps.(afield) = cat(2,unit_maps.(afield),unit_maps_sh{n}.(afield));
    end
end
%%% inclusion criteria
if strcmp(unittype,'all')
    includeflag = unit_maps.have_retin_flag ;
elseif strcmp(unittype,'KSgood')
    includeflag = unit_maps.goodKS_flag & unit_maps.have_retin_flag;
elseif strcmp(unittype,'taskmatched')
    includeflag = unit_maps.matchedunit_flag & unit_maps.have_retin_flag;
end


% 3Dmaps are raw response or predicted from peak, in XY grid, third
% dimension being the unit number. If Normalize, each unit is divided by
% its peak.
Pmap_ThreeDMat = cat(3,unit_maps.predictedmap{includeflag});
Pmap_Peaks = repmat(max(Pmap_ThreeDMat,[],[1,2]),[size(Pmap_ThreeDMat,1),size(Pmap_ThreeDMat,2)]);
Respmap_ThreeDMat = cat(3,unit_maps.RespMag{includeflag});
Resmap_Peaks = repmat(max(Respmap_ThreeDMat,[],[1,2]),[size(Respmap_ThreeDMat,1),size(Respmap_ThreeDMat,2)]);
if Normalize
    Pmap_ThreeDMat = Pmap_ThreeDMat./Pmap_Peaks;
    Respmap_ThreeDMat = Respmap_ThreeDMat./Resmap_Peaks;
end

%%% average: weighted or unweighted
imagesc(h1,nanmean(Respmap_ThreeDMat,3))
h1.Title.String =[areaname,'-', unittype];
im=imagesc(h2,[0.5,5.5],[0.5,4.5],nanmean(Pmap_ThreeDMat,3));
cmap = gray(10);% 
colormap(flipud(cmap)) %
im.AlphaData = 0.7; % 
axis tight
for cen = find(includeflag)
    hold on; scatter(unit_maps.xgcenter{(cen)},unit_maps.ygcenter{(cen)},36,'k.')% 10 or 36
end
h2.Title.String = [areaname,'-',unittype];


