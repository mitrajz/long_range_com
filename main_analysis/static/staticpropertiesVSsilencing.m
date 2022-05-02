% Before running this script, make sure retin_pipeline.m has been run for
% the animals, depth and spikeshape has been assiged, and
% makecells_generic* has been run for all included animals with retin,
% depth, and spikeshape assignment
%% IMPORTANT:to run this, load cells (def 150ms),  assign spikeshapes, set targetcell below, then run
% clean and remove according to the criteria: default: R0C1

targetcell =LMcells;

% smb or smb_ls
allgo = cell2mat(cellfun(@(x) x.smb_ls.go, targetcell,'UniformOutput',0));
allnogo = cell2mat(cellfun(@(x) x.smb_ls.nogo, targetcell,'UniformOutput',0));
% 
% allgo_s = allgo; allnogo_s = allnogo;
% allgo_s(abs(allgo_s)<gomeanconfs) = nan;
% allnogo_s(abs(allnogo_s)<nogomeanconfs) = nan;
% allgo = allgo_s; allnogo = allnogo_s;


allgobs = cell2mat(cellfun(@(x) x.smb_bs.go, targetcell,'UniformOutput',0));
allnogobs = cell2mat(cellfun(@(x) x.smb_bs.nogo, targetcell,'UniformOutput',0));

allretinx=cell2mat(cellfun(@(y) repmat(y,1,8),cellfun(@(x) {x.retin.center(1)}, targetcell),'UniformOutput',0));
allretiny=cell2mat(cellfun(@(y) repmat(y,1,8),cellfun(@(x) {x.retin.center(2)}, targetcell),'UniformOutput',0));
% 2 options: can change change from .distancelocal (distance to center of local shanks) to .distancelong 
% (distance to center of shanks in the other area)
allretinlocald=cell2mat(cellfun(@(y) repmat(y,1,8),cellfun(@(x) {x.retin.distancelong}, targetcell),'UniformOutput',0));

allAvfr=cell2mat(cellfun(@(y) repmat(y,1,8),cellfun(@(x) {x.Avfr},targetcell),'UniformOutput',0));
alldepth=cell2mat(cellfun(@(y) repmat(y,1,8),cellfun(@(x) {x.depth},targetcell),'UniformOutput',0));
allspikeamp=cell2mat(cellfun(@(y) repmat(y,1,8),cellfun(@(x) {x.spikeamplitude},targetcell),'UniformOutput',0));
allspikewidth=cell2mat(cellfun(@(y) repmat(y,1,8),cellfun(@(x) {x.spikewidth},targetcell),'UniformOutput',0));
% index of CENTERS
alltimecat = cell2mat(cellfun(@(x) {[1:8]},targetcell));


monitormiddlex = 3;
monitormiddley = 2.5;
allcenterdistance = sqrt( (allretinx - monitormiddlex).^2 + (allretiny - monitormiddley).^2 );

% figure;scatter([allcenterdistance,allcenterdistance],[allgo,allnogo],100,'k.')
% figure;scatter([allcenterdistance,allcenterdistance],[allgobs,allnogobs],100,'k.')
%% remove all nans and infs

nonanind = find(~isnan(mean([allretinlocald;allcenterdistance;allAvfr;alldepth;allspikeamp;allspikewidth],1)));
allretinlocald = allretinlocald(nonanind);
allcenterdistance = allcenterdistance(nonanind);
allAvfr = allAvfr(nonanind);
alldepth = alldepth(nonanind);
allspikeamp = allspikeamp(nonanind);
allspikewidth = allspikewidth(nonanind);
alltimecat = alltimecat(nonanind);


allgo = allgo(nonanind);
allnogo = allnogo(nonanind);
allgobs = allgobs(nonanind);
allnogobs = allnogobs(nonanind);

%% Example exploratory plots, relation of receptive field location and silencing

% figure;scatter(allretinx,allretiny,100,allnogo,'.');xlim([1,5]);ylim([1,4])
% figure;scatter([allretinlocald,allretinlocald],[allgo,allnogo],100,'k.')
% figure;subplot(1,2,1);scatter([allretinx,allretinx],[allgo,allnogo],100,'k.')
% subplot(1,2,2);scatter([allgo,allnogo],[allretiny,allretiny],100,'k.')
%%

rngnum = 1; % 2 or 1(def)
nfold = 20; %20
norm2obv = 1; % numobservations would be different when using different predictors, due to different nans
% Keep 1 in order to be able to compare.

glmprop.link = 'identity';
glmprop.spec = 'quadratic'; 
glmprop.cat = []; % index of categorical predictors: default: []

figure;

plotprop.plotprediction = 0;
plotprop.plotLL = 1;
plotprop.LLhandle = gca;
plotprop.LLloc = 1;
plotprop.singlectrl = 0; 
plotprop.scatterpl = 0;
plotprop.comparisons = 'perfold'; % 'perfold' or 'perobsv': discribed in methods. Perobs: not sure if complemety correct

tasktype = 'all'; % all, go, nogo


%%% prepare matrices:
clear X Y
Y.go = allgo; Y.nogo = allnogo;

X.go = allgobs; X.nogo = allnogobs;
plotprop.LLloc = 1;
[outmodels,Normalized_logl1] = predict_Slc_static(X,Y,tasktype,rngnum,nfold,norm2obv,glmprop,plotprop);
%%%%%%%%%%%%%%%% retin local distance -- or long range distance (depending on the first section)
X.go = [allgobs;allretinlocald]; X.nogo = [allnogobs;allretinlocald];
plotprop.LLloc = 2;
[outmodels,Normalized_logl2] = predict_Slc_static(X,Y,tasktype,rngnum,nfold,norm2obv,glmprop,plotprop);

if plotprop.singlectrl
    X.go = allretinlocald; X.nogo = allretinlocald;
    plotprop.LLloc =2.2;
    [outmodels,Normalized_logl] = predict_Slc_static(X,Y,tasktype,rngnum,nfold,norm2obv,glmprop,plotprop);
end
%%%%%%%%%%%%%%% retin center monitor distance
X.go = [allgobs;allcenterdistance]; X.nogo = [allnogobs;allcenterdistance];
plotprop.LLloc = 3;
[outmodels_r,Normalized_logl3] = predict_Slc_static(X,Y,tasktype,rngnum,nfold,norm2obv,glmprop,plotprop);

if plotprop.singlectrl
    X.go = allcenterdistance; X.nogo = allcenterdistance;
    plotprop.LLloc =3.2;
    [outmodels,Normalized_logl] = predict_Slc_static(X,Y,tasktype,rngnum,nfold,norm2obv,glmprop,plotprop);

end
%%%%%%%%%%%%%%% Avfr
X.go = [allgobs;allAvfr]; X.nogo = [allnogobs;allAvfr];
plotprop.LLloc = 4;
[outmodels,Normalized_logl4] = predict_Slc_static(X,Y,tasktype,rngnum,nfold,norm2obv,glmprop,plotprop);

if plotprop.singlectrl
    X.go = allAvfr; X.nogo = allAvfr;
    plotprop.LLloc = 4.2;
    [outmodels,Normalized_logl] = predict_Slc_static(X,Y,tasktype,rngnum,nfold,norm2obv,glmprop,plotprop);
end
%%%%%%%%%%%%%%% depth
plotprop.plotprediction = 0;

X.go = [allgobs;alldepth]; X.nogo = [allnogobs;alldepth];
plotprop.LLloc = 5;
[outmodels_d,Normalized_logl5] = predict_Slc_static(X,Y,tasktype,rngnum,nfold,norm2obv,glmprop,plotprop);

if plotprop.singlectrl
    X.go = alldepth; X.nogo = alldepth;
    plotprop.LLloc =5.2;
    [outmodels,Normalized_logl] = predict_Slc_static(X,Y,tasktype,rngnum,nfold,norm2obv,glmprop,plotprop);

end
plotprop.plotprediction = 0;


%%%%%%%%%%%%%%% spike width
X.go = [allgobs;allspikewidth]; X.nogo = [allnogobs;allspikewidth];
plotprop.LLloc = 6;
[outmodels,Normalized_logl6] = predict_Slc_static(X,Y,tasktype,rngnum,nfold,norm2obv,glmprop,plotprop);

if plotprop.singlectrl
    X.go = allspikewidth; X.nogo = allspikewidth;
    plotprop.LLloc =6.2;
    [outmodels,Normalized_logl] = predict_Slc_static(X,Y,tasktype,rngnum,nfold,norm2obv,glmprop,plotprop);

end
%%%%%%%%%%%%%%% spike amplitude
X.go = [allgobs;allspikeamp]; X.nogo = [allnogobs;allspikeamp];
plotprop.LLloc = 7;
[outmodels,Normalized_logl] = predict_Slc_static(X,Y,tasktype,rngnum,nfold,norm2obv,glmprop,plotprop);

if plotprop.singlectrl
    X.go = allspikeamp; X.nogo = allspikeamp;
    plotprop.LLloc =7.2;
     [outmodels,Normalized_logl] = predict_Slc_static(X,Y,tasktype,rngnum,nfold,norm2obv,glmprop,plotprop);

end


%%%%%%%%%%%%%%% categoical: tasktype: only when all
glmprop.cat = [];
X.go = [allgobs;zeros(size(allgobs))];X.nogo = [allnogobs;ones(size(allnogobs))];
plotprop.LLloc = 8;
[outmodels,Normalized_logl8] = predict_Slc_static(X,Y,tasktype,rngnum,nfold,norm2obv,glmprop,plotprop);
glmprop.cat = [];

if plotprop.singlectrl
    glmprop.cat = 1;
    X.go = zeros(size(allgobs)); X.nogo = zeros(size(allnogobs));
    plotprop.LLloc = 8.2;
    [poutmodels,Normalized_logl] = predict_Slc_static(X,Y,tasktype,rngnum,nfold,norm2obv,glmprop,plotprop);
    glmprop.cat = [];
end


%%%%%%%%%%%%%% catergorical: time
plotprop.plotprediction = 0;
glmprop.cat = [];
X.go = [allgobs;alltimecat]; X.nogo = [allnogobs;alltimecat];
plotprop.LLloc = 9;
[outmodels,Normalized_logl9] = predict_Slc_static(X,Y,tasktype,rngnum,nfold,norm2obv,glmprop,plotprop);
%glmprop.link = 'identity';
glmprop.cat = [];

if plotprop.singlectrl
    glmprop.cat = [];
    X.go = [alltimecat]; X.nogo = [alltimecat];
    plotprop.LLloc = 9.2;
    [outmodels,Normalized_logl] = predict_Slc_static(X,Y,tasktype,rngnum,nfold,norm2obv,glmprop,plotprop);
    glmprop.cat = [];
end
plotprop.plotprediction = 0;




%%%%%%%%%%%
ylim(plotprop.LLhandle,[-0.2 0.5])
xlim(plotprop.LLhandle,[0 plotprop.LLloc+1])
% set(plotprop.LLhandle,'XTick',[1 2 2.2 3 3.2 4 4.2 5 5.2 6 6.2 7 7.2 8 8.2 9 9.2])
% set(plotprop.LLhandle,'XTickLabel',{'bs','localDist+bs','localDist','centerDist+bs','centerDist',...
%     'Avfr+bs','Avfr','depth+bs','depth',...
%    'spikeWidth+bs','spikeWidth', 'spikeAmp+bs','spikeAmp',...
%     'task+bs','task','time+bs','time'})
% set(plotprop.LLhandle,'XTickLabelRotation',90)

set(plotprop.LLhandle,'XTick',[1 2 2.1 3 3.1 4 4.1 5 5.1 6 6.1 7 7.1 8 8.1 9 9.1])
set(plotprop.LLhandle,'XTickLabel',{'bs','longDist+bs','chance','centerDist+bs','chance',...
    'Avfr+bs','chance','depth+bs','chance',...
   'spikeWidth+bs','chance', 'spikeAmp+bs','chance',...
    'task+bs','chance','time+bs','chance'})
set(plotprop.LLhandle,'XTickLabelRotation',90)

%%%%


%% for a given 2 predictor model, see the shape of dependence predictors

allcoefs = cell2mat(cellfun(@(x) x.Coefficients.Estimate, outmodels_d,'UniformOutput',0))
%t=0:0.01:3; % set the range appropriately
t=-100:10:600; % set the range appropriately

% decide log or not!
tout = mean(allcoefs(1,:)) + t*(mean(allcoefs(3,:))) + (t.^2) *(mean(allcoefs(6,:)));
figure;plot(t,(tout))
