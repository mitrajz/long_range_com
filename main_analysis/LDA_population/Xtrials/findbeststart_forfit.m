function [bestmdl,beststart] = findbeststart_forfit(x,y,fo,ftype)

nrep = 200;%100;%1000, later 100
% sse is sum of square errors
meansse = nan(nrep,1);
spts = nan(nrep,3);
allmdls = cell(1,nrep);

% cross validation: (limit:number of folds: if using 80% trials, might fall under min number of trials)
% fixed seed 
numfolds = 5; 
rng(1,'twister');
cvind=crossvalind('Kfold',numel(x),numfolds);

%temp = [];
rng(1,'twister');
for i=1:nrep
    % A between 0 and 1
    % B between 0 and 1
    % tau between 0 and 1000
    
    fo.StartPoint= [rand rand 1000*rand];
    %fo.StartPoint= [.5 0.5 10*(i-1)+1];
    
    foldmeansse = nan(1,numfolds);
    for fold=1:numfolds
        traininds = (find(cvind ~= fold));
        testinds = (find(cvind == fold));
        
        ftype = setoptions(ftype,fo);
        try
            [mdl,gof,~] = fit(x(traininds),y(traininds),ftype);       
            yhat=feval(mdl,x(testinds));
            foldmeansse(fold)=sum((yhat-y(testinds)).^2)/numel(testinds);
        catch
           foldmeansse(fold) = nan;
        end
    end
    meansse(i) = mean(foldmeansse);
    allmdls{i} = mdl; 
    
    spts(i,:) = fo.StartPoint;

end

[~,ind] = min(meansse);
beststart = spts(ind,:);
bestmdl = allmdls{ind};