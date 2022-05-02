% checking th2, th1, prctmax: to see if spike shapes are classified nicely
% with these params
% just visualization, V1cells, LMcells dont change
%% params

%%
targetcell = eval(targetcellname);

clear spkind1 spkind2 

allspikes = cell2mat(cellfun(@(x) x.spiketemplate,targetcell,'UniformOutput',0));
allamplitudes = cell2mat(cellfun(@(x) x.spikeamplitude,targetcell,'UniformOutput',0));
allwidths = cell2mat(cellfun(@(x) x.spikewidth,targetcell,'UniformOutput',0));

spkind1 = find((allwidths)>th & (allwidths)<th2);
spkind2 = find((allwidths)<=th);

figure;scatter(allamplitudes,allwidths,'k')
hold on;scatter(allamplitudes(spkind1),allwidths(spkind1),'r')
hold on;scatter(allamplitudes(spkind2),allwidths(spkind2),'b')


figure;plot(allspikes(:,spkind1),'r')
hold on;plot(1+allspikes(:,spkind2),'b')
hold on;plot(2+allspikes,'k')


% spikeidth(still param dependent) : cell2mat(cellfun(@(x) (find(x.spiketemplate<min(prctmax*x.spiketemplate),1,'last') - find(x.spiketemplate<min(prctmax*x.spiketemplate),1,'first')),V1cells,'UniformOutput',0))