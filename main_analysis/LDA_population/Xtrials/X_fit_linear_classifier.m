function [vec,cnst,vec_n,mu_dist,mu_dist_norm,lda_norm_dist,Batta_dist] = X_fit_linear_classifier(lmodel,X,Y)



% X: nobsv * npredictors(cells)
% Y: nobsv * 1
% alltraces is Ntrials*Nneurons for both laser and baseline
Xi = X; Yi=Y;

vec = zeros(size(Xi,2),1); 
cnst = nan;
mu_dist = zeros(size(Xi,2),1); 
mu_dist_norm = nan;
vec_n = zeros(size(Xi,2),1);
lda_norm_dist = nan;
Batta_dist = nan;
%%% if shuffle: instead of laser trial: a bootstrapped sample of baseline
%%% trilas
if lmodel.shuffle == 1
    
    bsind = find(Yi==0);
    lsind = find(Yi==1);
    
     
    bsrand = bsind(randi(length(bsind),length(bsind),1));
    lsrand = lsind(randi(length(lsind),length(lsind),1));
    
      
    X = Xi([bsind;bsrand],:);
    Y = [Yi(bsind);Yi(bsrand)+1];
    
end



zerovars = union(find(isnan(mean(X,1))),find(abs(mean(X,1))<=lmodel.Xth));
nonzerovars = setdiff(1:size(X,2),zerovars);

% removing predictors (not observations)
X(:,zerovars)=[];
X = X - repmat(nanmean(X,1),[size(X,1),1]);


if strcmp(lmodel.type,'lda')
    try
        Mdlinear = fitcdiscr(X,Y,'Delta',lmodel.hyperparams.delta,'Gamma',lmodel.hyperparams.gamma,...
            'OptimizeHyperparameters','none','DiscrimType','Linear','Prior','uniform');% .2,.5
        vec(nonzerovars) =  Mdlinear.Coeffs(2,1).Linear;
        cnst = Mdlinear.Coeffs(2,1).Const;
        vec_n(nonzerovars) = vec(nonzerovars)/norm(vec(nonzerovars));
        
        % norm of the ldas vector. is relation to cov and deta mu: it is
        % equal to sqrt(Mudif*pinv(Mdlinear.Sigma^2)*Mudif') - if delta =
        % 0. The former coud be cauculated if we don't want delta to change
        % our result (gives similar result)
        lda_norm_dist = norm(Mdlinear.Coeffs(2,1).Linear); 
        
        %%% Battacharyya distance: this is the appropriate distance
        %%% between two distributions. This is based on the regularized
        %%% covariance (so gamma is already taken into account). sigma
        %%% is inverted using psuedoinverse. This based of the
        %%% assumtion that the 2 distribution have the same covariance
        %%% (pooled), same assumtion as LDA, and therefore the second
        %%% term is zero
        Mudif = Mdlinear.Mu(2,:) - Mdlinear.Mu(1,:);
        CovMat = Mdlinear.Sigma;
        Batta_dist = (Mudif*pinv(CovMat)*Mudif')/8;

    catch
        disp('caught exception')
    end
    
    % distance between class means
    try
        mu_dist(nonzerovars) = (Mdlinear.Mu(1,:)' - Mdlinear.Mu(2,:)');
        mu_dist_norm = norm(mu_dist); 
        mu_dist(nonzerovars) = mu_dist(nonzerovars)/norm(mu_dist(nonzerovars));
        
    catch
    end
    
    
    
    
elseif strcmp(lmodel.type,'svm')
    disp('svm not implemented yet')
end











