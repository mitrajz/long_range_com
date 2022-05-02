function [out_vec,out_train_er,out_test_er,proj,test_er,train_er] = X_lda_withCV(X,Y,lmodel,nfold,alltraces,nonzerovars,projvec)
LOOCV = 0;
usematlabpredict = 0;
% fix seed
rng(1,'twister');




Xraw = X;
Yraw = Y;
bsind = find(Yraw==0);
lsind = find(Yraw==1);


% making cv indices
cvind=crossvalind('Kfold',numel(bsind)+numel(lsind),nfold);

if LOOCV
    cvind_bs = bsind;
    cvind_ls = lsind;
    nfold = numel(cvind_bs) + numel(cvind_ls);
end

% initialization: return these
vec_n = zeros(size(alltraces,2),nfold);
vec = zeros(size(alltraces,2),nfold);
cnst = nan(1,nfold);
train_er = nan(1,nfold);
test_er = nan(1,nfold); 
% ntrials*time
proj = nan(size(alltraces,1),size(alltraces,3));
for fold=1:nfold
    traininds = (find(cvind ~= fold));
    testinds = (find(cvind == fold));
            
    X = Xraw(traininds,:);
    Y = Yraw(traininds);
    Xtest = Xraw(testinds,:);
    Ytest = Yraw(testinds);
    try
         Mdlinear = fitcdiscr(X,Y,'Delta',lmodel.hyperparams.delta,'Gamma',lmodel.hyperparams.gamma,...
             'OptimizeHyperparameters','none','DiscrimType','Linear','Prior','uniform');% .2,.5

       vec(nonzerovars,fold) =  Mdlinear.Coeffs(2,1).Linear;
       cnst(fold) = Mdlinear.Coeffs(2,1).Const;
        vec_n(nonzerovars,fold) = vec(nonzerovars,fold)/norm(vec(nonzerovars,fold));
        %%% claculate errors
        %%% onlt nonzerovars cells are used to calculate error
       
        if ~usematlabpredict
            % train
            Pred = [ones(size(X,1),1),X] * [cnst(fold);vec(nonzerovars,fold)];
            Pred = Pred>0;
            train_er(fold) = numel(find(xor(Pred,Y)))/numel(Y);
            % test
            Pred = [ones(size(Xtest,1),1),Xtest] * [cnst(fold);vec(nonzerovars,fold)];
            Pred = Pred>0;
            test_er(fold) = numel(find(xor(Pred,Ytest)))/numel(Ytest);
        else
            %Pred = predict(Mdlinear,Xtest);
        end
        
       if lmodel.doproj == 1
           for tr = 1:length(testinds)
             proj(testinds(tr),:) = (squeeze(alltraces(testinds(tr),:,:))'*vec_n(:,fold))';
           end
       elseif lmodel.doproj == 2
           for tr = 1:length(testinds)
             proj(testinds(tr),:) = (squeeze(alltraces(testinds(tr),:,:))'*projvec)';
           end
       end
    catch
        disp('caught exception')
    end
end

%

out_vec = nanmean(vec,2);
out_train_er = nanmean(train_er);
out_test_er = nanmean(test_er);