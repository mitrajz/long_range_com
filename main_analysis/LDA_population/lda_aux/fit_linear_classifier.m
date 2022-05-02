function [vec,cnst,vec_n,mu_dist,cv,btstrp,proj,test_er,train_er] = fit_linear_classifier(X,Y,lmodel,alltraces,projvec)


% X: nobsv * npredictors(cells)
% Y: nobsv * 1
% alltraces is Ntrials*Nneurons for both laser and baseline

% params
nboot =100;
if isfield(lmodel,'cvn')
    nfold = lmodel.cvn;
else
    nfold = 10;
end
% initialization
cv = struct;
btstrp = struct;
proj = nan(size(alltraces,1),size(alltraces,3));
Xi = X; Yi=Y;

vec = zeros(size(alltraces,2),1);
cnst = nan;
mu_dist = nan(size(alltraces,2),1); 
vec_n = zeros(size(alltraces,2),1);

% 
zerovars = union(find(isnan(mean(X,1))),find(abs(mean(X,1))<=lmodel.Xth));
nonzerovars = setdiff(1:size(X,2),zerovars);

% removing predictors (not observations)
X(:,zerovars)=[];
X = X - repmat(nanmean(X,1),[size(X,1),1]);


if strcmp(lmodel.type,'lda')
    if ~lmodel.cv
        Mdlinear = fitcdiscr(X,Y,'Delta',lmodel.hyperparams.delta,'Gamma',lmodel.hyperparams.gamma,...
            'OptimizeHyperparameters','none','DiscrimType','Linear','Prior','uniform');% .2,.5
        %ldacoeffs_go = [Mdlinear.Coeffs(2,1).Const ; Mdlinear.Coeffs(2,1).Linear];
        vec(nonzerovars) =  Mdlinear.Coeffs(2,1).Linear;
        cnst = Mdlinear.Coeffs(2,1).Const;
        %%%%% **** THIS PROJ IS ON VEC NOT VEC_N
        if lmodel.doproj == 1
            for tr = 1:size(alltraces,1)
                proj(tr,:) = (squeeze(alltraces(tr,:,:))'*vec)';
            end
        elseif lmodel.doproj ==2
            for tr = 1:size(alltraces,1)
                proj(tr,:) = (squeeze(alltraces(tr,:,:))'*projvec)';
            end
        end
    else
        [vec,cv.train_er,cv.test_er,proj,test_er,train_er] = lda_withCV(X,Y,lmodel,nfold,alltraces,nonzerovars,projvec);
    end
    
    %
    % distance between class means
    try
        mu_dist(nonzerovars) = (Mdlinear.Mu(1,:)' - Mdlinear.Mu(2,:)');
    catch
    end
elseif strcmp(lmodel.type,'svm')
    CVMdl = fitclinear(X,Y,'ObservationsIn','rows',...
        'Learner','svm','Regularization','ridge',...
        'Lambda',lmodel.hyperparams.gamma,'GradientTolerance',1e-8);
    vec(nonzerovars) = CVMdl.Beta;
end
% normalized vector
vec_n = vec/norm(vec);
%% bootstrap
if lmodel.boot
    %%% If bootstrap, same code as above, but with resampled trial indices.
    %%% Only the bootstrapped  vectors are saves (in btstrp fields), and CV
    %%% errors, proj,.. are discarded. Bootstrapping happens with the same
    %%% model specification as above, eg. if lmode.CV == 1, bootstrap also
    %%% resamples the crossvalidationn model 
    
    Xraw = Xi;
    Yraw = Yi;
    bsind = find(Yi==0);
    lsind = find(Yi==1);
    % initialization
    btstrp.vec = zeros(size(Xraw,2),nboot); 
   % btstrp.cnst = nan(1,nboot);
    btstrp.mu_dist = zeros(size(Xraw,2),nboot); 
    btstrp.vec_n = zeros(size(Xraw,2),nboot);
    btstrp.cv = struct;
    btstrp.cv.train_er = cell(1,0);
    btstrp.cv.test_er = cell(1,0);
    btstrp.X = cell(1,nboot);
    btstrp.Y = cell(1,nboot);
    
    for i=1:nboot
        rng(i,'twister');
        bsrand = bsind(randi(length(bsind),length(bsind),1));
        lsrand = lsind(randi(length(lsind),length(lsind),1));
        if lmodel.shuffle == 0
            X = Xraw([bsrand;lsrand],:);
            Y = Yraw([bsrand;lsrand]);
        elseif lmodel.shuffle == 1
             X = Xraw;
             Y = Yraw(randperm(length(Yraw)));
        elseif lmodel.shuffle == 2
            X = Xraw([bsind;bsrand],:);
            Y = [Yraw(bsind);Yraw(bsrand)+1];
        end  
        btstrp.X{i} = X;
        btstrp.Y{i} = Y;
       
        % 
        zerovars = union(find(isnan(mean(X,1))),find(abs(mean(X,1))<=lmodel.Xth));
        nonzerovars = setdiff(1:size(X,2),zerovars);
        % removing predictors (not observations)
        X(:,zerovars)=[];
        X = X - repmat(mean(X,1),[size(X,1),1]);
        if strcmp(lmodel.type,'lda')
            if ~lmodel.cv
                Mdlinear = fitcdiscr(X,Y,'Delta',lmodel.hyperparams.delta,'Gamma',lmodel.hyperparams.gamma,...
                    'OptimizeHyperparameters','none','DiscrimType','Linear','Prior','uniform');% .2,.5
                %ldacoeffs_go = [Mdlinear.Coeffs(2,1).Const ; Mdlinear.Coeffs(2,1).Linear];
                btstrp.vec(nonzerovars,i) =  Mdlinear.Coeffs(2,1).Linear;
              %  btstrp.cnst(i) = Mdlinear.Coeffs(2,1).Const;
            else
                [ btstrp.vec(:,i),~,~,~,btstrp.cv.test_er{end+1},btstrp.cv.train_er{end+1}] = ...
                    lda_withCV(X,Y,lmodel,nfold,alltraces,nonzerovars,projvec);
            end
            
            % distance between class means
            try
                btstrp.mu_dist(nonzerovars,i) = (Mdlinear.Mu(1,:)' - Mdlinear.Mu(2,:)');
            catch
            end
            
        elseif strcmp(lmodel.type,'svm')
            CVMdl = fitclinear(X,Y,'ObservationsIn','rows',...
                'Learner','svm','Regularization','ridge',...
                'Lambda',lmodel.hyperparams.gamma,'GradientTolerance',1e-8);
            btstrp.vec(nonzerovars,i) = CVMdl.Beta;
        end
        % normalized vector
        btstrp.vec_n(:,i) = btstrp.vec(:,i)/norm(btstrp.vec(:,i));
    end
    
end

