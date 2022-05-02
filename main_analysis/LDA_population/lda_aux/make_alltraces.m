function proj = make_alltraces(X,targetcell,lag,text,lmodel)
if  lmodel.eqtrials
    if strcmp(text,'go')
        % alltrials*ncells*time
        proj = nan(size(X,1),size(X,2),size(targetcell{1}.nbs.go{lag},2));
        for j = 1:size(X,2)
            % nalltrials*timr
            proj(:,j,:)=[targetcell{j}.nbs.go{lag};targetcell{j}.nls.go{lag}];
        end
        
    elseif strcmp(text,'nogo')
        % alltrials*ncells*time
        proj = nan(size(X,1),size(X,2),size(targetcell{1}.nbs.nogo{lag},2));
        for j = 1:size(X,2)
            % nalltrials*timr
            proj(:,j,:)=[targetcell{j}.nbs.nogo{lag};targetcell{j}.nls.nogo{lag}];
        end
        
    else
        error('something wrong')
    end
else
    proj = nan(size(X,1),size(X,2),size(targetcell{1}.nbs.go,2));
end