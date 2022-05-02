%% cross validated ordering

targetcell = V1cells;
alleffects_train = nan(length(targetcell),nlags);
alleffects_test = nan(length(targetcell),nlags);
rng(1) % fix at any number
for cellnum = 1:numel(targetcell)   
    
   nummintrials = min([cellfun(@(x) numel(x), targetcell{cellnum}.laAls.go),...
            cellfun(@(x) numel(x), targetcell{cellnum}.laAbs.go)]);
        trialrandord = randperm(nummintrials);
        trainind = trialrandord(1:round(.8*nummintrials));
        testind = trialrandord(round(.8*nummintrials)+1:nummintrials);
        for l = 1:nlags
           alleffects_train(cellnum,l) = ...
                (cell2mat(cellfun(@(x) nanmean(x(trainind)),targetcell{cellnum}.laAls.go(l),'UniformOutput',0)) - ...
                 cell2mat(cellfun(@(x) nanmean(x),targetcell{cellnum}.laAbs.go(l),'UniformOutput',0)))/...
                 cell2mat(cellfun(@(x) nanmean(x),targetcell{cellnum}.laAbs.go(l),'UniformOutput',0));
            
               alleffects_test(cellnum,l) = ...
               (cell2mat(cellfun(@(x) nanmean(x(testind)),targetcell{cellnum}.laAls.go(l),'UniformOutput',0)) - ...
                 cell2mat(cellfun(@(x) nanmean(x),targetcell{cellnum}.laAbs.go(l),'UniformOutput',0)))/...
                 cell2mat(cellfun(@(x) nanmean(x),targetcell{cellnum}.laAbs.go(l),'UniformOutput',0));
                      
        end
end
%%% 
alleffects_train(find(isinf(alleffects_train))) = 0;
alleffects_train(find(isnan(alleffects_train))) = 0;
alleffects_train(abs(cell2mat(cellfun(@(x) x.smb_z.go,targetcell,'uniformoutput',0)'))<2)=0;

alleffects_test(find(isinf(alleffects_test))) = 0;
alleffects_test(find(isnan(alleffects_test))) = 0;
alleffects_test(abs(cell2mat(cellfun(@(x) x.smb_z.go,targetcell,'uniformoutput',0)'))<2)=0;


[~,ind]=max(alleffects_train');
[~,order]=sort(ind);
figure;
alleffects_train_sorted = alleffects_train(order,:);
subplot(1,2,1);imagesc(alleffects_train_sorted);set(gca,'Clim',[-1 1]);colormap('parula');
alleffects_test_sorted = alleffects_test(order,:);


alleffects_test_sorted((sum(alleffects_test_sorted,2) == 0),:) = [];
subplot(1,2,2);imagesc(alleffects_test_sorted);set(gca,'Clim',[-1 1]);colormap('parula');

