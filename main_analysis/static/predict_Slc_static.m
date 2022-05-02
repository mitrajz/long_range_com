function [outmodels,Normalized_logl] = predict_Slc_static(Xs,Ys,tasktype,rngnum,nfold,norm2obv,glmprop,plotprop)

%% GLM predincting silencing activity(smb_slc) based on retininfo

if strcmp(tasktype,'go')
    X = Xs.go;
    Y = Ys.go;
elseif strcmp(tasktype,'nogo')
    X = Xs.nogo;
    Y = Ys.nogo;
elseif strcmp(tasktype,'all')
    X = [Xs.go,Xs.nogo];
    Y = [Ys.go,Ys.nogo];
end
% removing imdexes that are nan in any of preditors or output(Y)
anynan = Y;
% for i=1:size(X,1)
%     anynan = anynan+X(i,:);
% end
nonnanindc = find(~isnan(anynan));
X = X(:,nonnanindc);
Y = Y(:,nonnanindc);
% removing imdexes that areinf in any of preditors or output(Y)
anyinf = Y;
% for i=1:size(X,1)
%      anyinf = anyinf+X(i,:);
% end
noninfindc = find(~isinf(anyinf));
X = X(:,noninfindc);
Y = Y(:,noninfindc);

% if model is log, log(y) = f(x), so Y can't be zero. for log: Y>0
if strcmp(glmprop.link,'log') || strcmp(glmprop.link,'logit') || ...
        strcmp(glmprop.link,'comploglog')
    X = X(:,find(Y~=0));
    Y = Y(:,find(Y~=0));
end
%%%
rng(rngnum)
cvind=crossvalind('Kfold',size(Y,2),nfold);



outmodels = cell(1,nfold);
Normalized_logl = [];
Normalized_logl_chance = [];



alltestX = [];
alltestY = [];
allYhat = [];

for fold=1:nfold
    traininds = (find(cvind ~= fold));
    testinds = (find(cvind == fold));
    
    %%%%%%%%%%%%% go
    % make main model and null model from training set
    mdl = fitglm(X(:,traininds)',Y(1,traininds)',glmprop.spec,'Link',glmprop.link,'CategoricalVars',glmprop.cat);
    outmodels{fold} = mdl;
    % log likelihood of test set under main model and test model
    testx = X(:,testinds)';
    testy = Y(1,testinds)';
    yhat = predict(mdl,testx);
    N = length(find(~isnan(yhat)));
    % yhatnull = predict(mdl,zeros(size(testx)));%repmat(mdl.Coefficients.Estimate(1),numel(yhat),1);%predict(mdl,zeros(size(testx)));%
    % sigmanull = std(yhatnull);
    sigma = sqrt(mdl.Dispersion); %std(yhat);% std(testy);
    

 
    %%%%%%%%%%%%% null model: only baseline trials: assuming first row is
    %%%%%%%%%%%%% baseline
    
    mdl_null = fitglm(X(1,traininds)',Y(1,traininds)',glmprop.spec,'Link',glmprop.link);
    testx_null = X(1,testinds(find(~isnan(yhat))))';
    testy_null = Y(1,testinds(find(~isnan(yhat))))';
    yhat_null = predict(mdl_null,testx_null);
    N_null = length(find(~isnan(yhat_null)));
    sigma_null = sqrt(mdl_null.Dispersion); %std(yhat);% std(testy);
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % saving real and predicted values, to check predictions.
    alltestX = [alltestX;testx];
    alltestY = [alltestY;testy];
    allYhat = [allYhat;yhat];
    
    %%%%%%%%% getting likelihood values
    [N N_null]
    if strcmp(plotprop.comparisons,'perfold') % over folds
         Normalized_logl =[Normalized_logl; ( ((-N/2)*Log2(2*pi*sigma^2)) - Log2(exp((sum((yhat - testy).^2)/(2*sigma^2))))   )/N - ...
            ( ((-N_null/2)*Log2(2*pi*sigma_null^2)) -  Log2(exp((sum((yhat_null - testy_null).^2)/(2*sigma_null^2))))   )/N_null];
    elseif strcmp(plotprop.comparisons,'perobsv') %  if per observation: chi^2
        if fold == 1 % over 1 test-train set
           Normalized_logl =[Normalized_logl; ((((-1/2)*Log2(2*pi*sigma^2)) - Log2(exp((((yhat - testy).^2)/(2*sigma^2))))  )) - ...
                (((-N_null/2)*Log2(2*pi*sigma_null^2)) - Log2(exp((sum((yhat_null - testy_null).^2)/(2*sigma_null^2))))   )/N_null];
            Normalized_logl_chance =[Normalized_logl_chance; (((-1/2)*Log2(2*pi*sigma_null^2)) - Log2(exp(( ((yhat_null - testy_null).^2)/(2*sigma_null^2))))  ) - ...
                (((-N_null/2)*Log2(2*pi*sigma_null^2)) - Log2(exp((sum((yhat_null - testy_null).^2)/(2*sigma_null^2))))  )/N_null];
        end
    end
end

% Normalized_logl is normalized to the number of observations.




%% plotting
% CV predictions vs values
if plotprop.plotprediction
    f=figure;
    for i=1:size(X,1)
        axp = subplot(1,size(X,1),i);%axes(f);
        scatter(axp,alltestX(:,i),alltestY,100,'k.');
        hold(axp,'on'); scatter(axp,alltestX(:,i),allYhat,100,'c.');
        axp.Title.String = (sprintf('predictor %d',i));
        
    end
end

% log likelihood plots
if plotprop.plotLL
    % boxplots. It doesnt support single variable model contrls at the
    % moment
    if plotprop.LLloc == round(plotprop.LLloc)
        if strcmp(plotprop.comparisons,'perfold')
            
            hold(plotprop.LLhandle,'on')
            if plotprop.scatterpl
                scatter(plotprop.LLhandle,zeros(size(Normalized_logl)) + plotprop.LLloc,Normalized_logl,'.','MarkerEdgeColor',[0 0 0],'MarkerEdgeAlpha',0.5)
                hold(plotprop.LLhandle,'on')
                scatter(plotprop.LLhandle,plotprop.LLloc,nanmean(Normalized_logl),'+','MarkerEdgeColor',[1 0 0],'MarkerEdgeAlpha',1)
                
            end
            
            hold(plotprop.LLhandle,'on')
            boxplot(plotprop.LLhandle,Normalized_logl,....
                zeros(size(Normalized_logl)) + plotprop.LLloc,'Positions',zeros(size(Normalized_logl)) + plotprop.LLloc,...
                'PlotStyle','compact','MedianStyle','line','OutlierSize',1,'Colors',[0 0 0]);
            hold(plotprop.LLhandle,'on');
            line([0 10],[0 0],'Color','k');
            hold(plotprop.LLhandle,'on');
            
            % the result of sign rank test would be exactly the same if instead of
            % normalizing to baseline per fold, calculate model and
            % nullmodel loglikelihoods and do a paired test between model
            % and nullmodel. Here, we normalize per fold and do a
            % signrank, comparison to zero med. Signrank instead of ttest
            % because distribution is not normal (chi-squared)
            text(plotprop.LLloc,0.4,num2str(signrank(Normalized_logl,0,'tail','right')),...
                'FontSize',12)
            
        elseif strcmp(plotprop.comparisons,'perobsv')
            hold(plotprop.LLhandle,'on')
            if plotprop.scatterpl
                scatter(plotprop.LLhandle,zeros(size(Normalized_logl)) + plotprop.LLloc,Normalized_logl,'.','MarkerEdgeColor',[0 0 0],'MarkerEdgeAlpha',0.5)
                hold(plotprop.LLhandle,'on')
                scatter(plotprop.LLhandle,plotprop.LLloc,nanmean(Normalized_logl),'+','MarkerEdgeColor',[1 0 0],'MarkerEdgeAlpha',1)
                
            end
            
            hold(plotprop.LLhandle,'on')
            boxplot(plotprop.LLhandle,Normalized_logl,....
                zeros(size(Normalized_logl)) + plotprop.LLloc,'Positions',zeros(size(Normalized_logl)) + plotprop.LLloc,...
                'PlotStyle','compact','MedianStyle','line','OutlierSize',1,'Colors',[0 0 0]);
            %%%%%% chance
            hold(plotprop.LLhandle,'on')
            if plotprop.scatterpl
                scatter(plotprop.LLhandle,0.1+zeros(size(Normalized_logl_chance)) + plotprop.LLloc,Normalized_logl_chance,'.','MarkerEdgeColor',[0.8 0.8 0.8],'MarkerEdgeAlpha',0.5)
                hold(plotprop.LLhandle,'on')
                scatter(plotprop.LLhandle,0.1+plotprop.LLloc,nanmean(Normalized_logl_chance),'+','MarkerEdgeColor',[1 0 0],'MarkerEdgeAlpha',1)
                
            end
            
            hold(plotprop.LLhandle,'on')
            boxplot(plotprop.LLhandle,Normalized_logl_chance,....
                0.1 + zeros(size(Normalized_logl_chance)) + plotprop.LLloc,'Positions',0.1 + zeros(size(Normalized_logl_chance)) + plotprop.LLloc,...
                'PlotStyle','compact','MedianStyle','line','OutlierSize',1,'Colors',[0.8 0.8 0.8]);
            
            hold(plotprop.LLhandle,'on');
            
            text(plotprop.LLloc,1,num2str(ranksum(Normalized_logl,Normalized_logl_chance)),...
                'FontSize',6)
        end
    end
    
    
    
end

    
end


