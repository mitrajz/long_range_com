function unit_maps = retin_process_unit(AllSpikes,clusters,Order,Window,PDAllOn,filt_pd,on_off,sres,ResponseOnset_ms)
% add: give reponse and baseline durations
% Here we do for all clusters, no checking of KSlabel or matching
% Response window is from ResponseOnset:end(200ms)
% baseline is pre zero, symmetric window

allclusterids = unique(clusters);

unit_maps = struct;
% cell with size of numberof units
unit_maps.clusterid = cell(1,length(allclusterids));
unit_maps.RespMag = cell(1,length(allclusterids));
unit_maps.xgcenter = cell(1,length(allclusterids));
unit_maps.ygcenter = cell(1,length(allclusterids));
unit_maps.p = cell(1,length(allclusterids));
unit_maps.predictedmap = cell(1,length(allclusterids));
unit_maps.predictedmapcoord = cell(1,length(allclusterids));




for cluster=1:length(allclusterids)
    
    unit_maps.clusterid{cluster} = allclusterids(cluster);
    
    Sp=nan(size(filt_pd));
    Sp(AllSpikes(find(clusters==cluster)))=1;
    
   
    
    tres=Window+1+ResponseOnset_ms*30;
    tbase=Window+1;
    
    % on 
    R0 = nan(2*Window+1,length(unique(Order)));
    RespMag_on = nan(1,4*5);
    for u= 1:(length(unique(Order))/2)
        nu=u;
        R0(1:2*Window+1,nu)= nansum(Sp(PDAllOn(find(Order==u),:))',2);
        RespMag_on(nu)= nansum(R0(tres:end,nu),1) - nansum(R0((tbase-length(R0(tres:end,nu))+1):tbase,nu),1);
    end
    RespMag_on = reshape(RespMag_on,5,4);
    RespMag_on = RespMag_on';
    % off
    R0 = nan(2*Window+1,length(unique(Order)));
    RespMag_off = nan(1,4*5);
    for u=(length(unique(Order))/2+1):length(unique(Order))
        nu=u-(length(unique(Order))/2);
        R0(1:2*Window+1,nu)= nansum(Sp(PDAllOn(find(Order==u),:))',2);
        RespMag_off(nu)= nansum(R0(tres:end,nu),1) - nansum(R0((tbase-length(R0(tres:end,nu))+1):tbase,nu),1);
    end
    
    RespMag_off = reshape(RespMag_off,5,4);
    RespMag_off = RespMag_off';
    %%% assigning Respmag field
    if isnan(on_off)
        unit_maps.RespMag{cluster} = (RespMag_on+RespMag_off)/2;
    elseif on_off
        unit_maps.RespMag{cluster} = RespMag_on;
    elseif ~on_off
        unit_maps.RespMag{cluster} = RespMag_off;
    end
    
    if max(max(unit_maps.RespMag{cluster})) >= 10 % if enough spikes in the maximum block
        px=nan(1,10);
        py=nan(1,10);
        value=nan(1,10);
        p=nan(10,6);
        for rep=1:1:10
            [px(rep), py(rep), value(rep),p(rep,1:6)] = fitgaussian2(unit_maps.RespMag{cluster});
        end
        
        [~,iv]=min(value);
        unit_maps.xgcenter{cluster} = px(iv);
        unit_maps.ygcenter{cluster} = py(iv);
        unit_maps.p{cluster} = p(iv,:);
                
        yvec = repmat(0.5:sres:4.5,1,length(0.5:sres:5.5))';
        xvec = repmat(0.5:sres:5.5,length(0.5:sres:4.5),1);
        xvec = xvec(:);
        
        unit_maps.predictedmapcoord{cluster}.x = xvec;
        unit_maps.predictedmapcoord{cluster}.y = yvec;
        
        
        predfuncplot = @(b,xvec,yvec) b(1) + b(2) .* exp(-(xvec - b(3)).^2./(b(4).^2)) .* exp(-(yvec - b(5)).^2./(b(6).^2));
        predicted = predfuncplot(unit_maps.p{cluster},unit_maps.predictedmapcoord{cluster}.x,unit_maps.predictedmapcoord{cluster}.y);
        unit_maps.predictedmap{cluster} = reshape(predicted,length(0.5:sres:4.5),length(0.5:sres:5.5));
            
        
    end
    
end
