function rep = cattorep(rep,trinds,i0,D,p1,p2,i,type)
param = 2; % subtract first half of stimulus

if strcmp (type,'go')
    for tr = trinds
        % one trials projection: length: number of time points
        % vec_n{i} or vec_n{4} determines auto or not projection
        
        rawroj = D.part{p1}.go.vec_n{i0}'*...
            squeeze(D.part{p2}.go.alltraces{i}(tr,:,:));
        
        % sbtract pre stimulus
        proj = (rawroj-nanmean(rawroj(:,1:size(rawroj,2)/param),2));

        rep(end+1,:) = proj;
        
    end
elseif strcmp (type,'nogo')
    for tr = trinds
        % one trials projection: length: number of time points
        % vec_n{i} or vec_n{4} determines auto or not projection
        
        rawroj = D.part{p1}.nogo.vec_n{i0}'*...
            squeeze(D.part{p2}.nogo.alltraces{i}(tr,:,:));
        
        % sbtract pre stimulus
        proj = (rawroj-nanmean(rawroj(:,1:size(rawroj,2)/param),2));

        rep(end+1,:) = proj;
        
    end
else
    error('check here')
end



