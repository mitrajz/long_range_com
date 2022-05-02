function [proj_vec_go,proj_vec_nogo] = get_proj_vec(targetcell,lmodel,i)
proj_vec_go = nan;
proj_vec_nogo = nan;
if lmodel.doproj == 2   
    lmodel.doproj = 1; 
    X = [cell2mat(cellfun(@(x) x.laAbs.go{i}',targetcell,'UniformOutput',0)')';...
        cell2mat(cellfun(@(x) x.laAls.go{i}',targetcell,'UniformOutput',0)')'];
    Y = [zeros(size(targetcell{i}.laAbs.go{i},1),1);ones(size(targetcell{i}.laAls.go{i},1),1)];
    % alltrials*ncells*time
    alltraces = make_alltraces(X,targetcell,i,'go',lmodel);
    % do not manipulate X,Y before this point (proj depends on it)
    proj_vec_go = zeros(1,size(X,2));
    [~,~,proj_vec_go,~,~,~] = fit_linear_classifier(X,Y,lmodel,alltraces,nan);
    
    X = [cell2mat(cellfun(@(x) x.laAbs.nogo{i}',targetcell,'UniformOutput',0)')';...
        cell2mat(cellfun(@(x) x.laAls.nogo{i}',targetcell,'UniformOutput',0)')'];
    Y = [zeros(size(targetcell{i}.laAbs.nogo{i},1),1);ones(size(targetcell{i}.laAls.nogo{i},1),1)];
    % alltrials*ncells*time
    alltraces = make_alltraces(X,targetcell,i,'nogo',lmodel);
    % do not manipulate X,Y before this point (proj depends on it)
    proj_vec_nogo = zeros(1,size(X,2));
    [~,~,proj_vec_nogo,~,~,~] = fit_linear_classifier(X,Y,lmodel,alltraces,nan);
end
