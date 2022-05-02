
logl = [];
for i=1:50

    targetcell = LMcells(i);
    
    allgo = cell2mat(targetcell{1}.laAls.go');
    allgobsav = [];
    alllag = [];
    alltime = [];
    for l = 1:8
        % version 1 : frm average basseline
        allgobsav = [allgobsav;repmat(nanmean(V1cells{i}.laAbs.go{l}),length(V1cells{i}.laAls.go{l}),1)];
        % version 2: from laser_ctrl (laser trials, one bin before laser)
        % allgobsav = [allgobsav;targetcell{1}.laAls_ctrl.go{l}];
        alllag = [alllag;repmat(l,length(targetcell{1}.laAls.go{l}),1)];
        alltime = [alltime;(1:length(targetcell{1}.laAls.go{l}))'];
    end
    
    rngnum = 2; 
    nfold = 20; 
    norm2obv = 1;
    
    glmprop.link = 'identity';
    glmprop.spec = 'purequadratic'; 
    glmprop.cat = []; 
    
   % figure;
    plotprop.plotprediction = 0;
    plotprop.plotLL = 0;
    plotprop.LLhandle = gca;
    plotprop.LLloc = 0;
    plotprop.singlectrl = 0; 
    plotprop.scatterpl = 0;
    plotprop.comparisons = 'perfold'; 
    tasktype = 'go';
    
    clear X Y
    Y.go = allgo';     
    X.go = [allgobsav';alllag'];
    glmprop.cat = [2];
    [outmodels,Normalized_logl1] = predict_Slc_static(X,Y,tasktype,rngnum,nfold,norm2obv,glmprop,plotprop);
    logl = [logl;Normalized_logl1];
    
    
end