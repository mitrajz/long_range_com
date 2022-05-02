% in hyper param sweep, projs and some other output are completely
% discarded. The only relevant variable pretty much is errors
% eqtrials is important param. Must be saved.
%% linear model properties
lmodel.type = 'lda';
lmodel.cv =1; % KEEP 1 for SWEEP.
lmodel.boot = 0; % KEEP 0 for SWEEP. if 1, model returns a nonempty btstrp struct with bootstrapped vectors
lmodel.Xth =0; % def : 0 
lmodel.eqtrials = 1;
lmodel.doproj = 0; % KEEP 0 for SWEEP 0: no proj,1: autoproj(to the trials not used to fit model), 2:to a given lag
lmodel.verbose = 0; %KEEP 0 for SWEEP

lmodel.hyperparams.gammarange = linspace(0,1,20);
lmodel.hyperparams.deltarange = linspace(0,1,20);

lmodel.exptype = 'FF';
lmodel.clean = 1;
lmodel.remove = 1;
lmodel.cvn = 10;
%%
% load data, clean and remove cells: (aggregate analysis)
% choose animal

numanimals = max([cellfun(@(x) x.simulcode, V1cells), cellfun(@(x) x.simulcode, LMcells)]);
nlags=8;
savepath = '/mnt/data/Mitra/figs/P2_L/bothdircombined/2020_02_24/LDA';
savetofile = 1;


%%% initialization
animalmodel = cell(1,numanimals);

vec_go = cell(1,nlags);
cnst_go = cell(1,nlags);
vec_n_go = cell(1,nlags);
cv_go = cell(1,nlags);

vec_nogo = cell(1,nlags);
cnst_nogo = cell(1,nlags);
vec_n_nogo = cell(1,nlags);
cv_nogo = cell(1,nlags);

proj_vec_go = nan;
proj_vec_nogo = nan;


for animalnum=1:numanimals
    if strcmp(lmodel.exptype,'FF')
        targetcell = LMcells(find(cellfun(@(x) x.simulcode,LMcells) == animalnum ));
    elseif strcmp(lmodel.exptype,'FB')
        targetcell = V1cells(find(cellfun(@(x) x.simulcode,V1cells) == animalnum ));
    end
      
    animalmodel{animalnum} = struct;
    animalmodel{animalnum}.animalname = unique(cellfun(@(x) x.animal, targetcell,'UniformOutput',0));
    animalmodel{animalnum}.MeanTestEr = nan(length(lmodel.hyperparams.gammarange),length(lmodel.hyperparams.deltarange));
    animalmodel{animalnum}.MeanTrainEr = nan(length(lmodel.hyperparams.gammarange),length(lmodel.hyperparams.deltarange));
    animalmodel{animalnum}.cleanTh = cleanupcriteria.smb_activity;
    %%%
    % equalize trials, if needed.
    if lmodel.eqtrials
        randeq = 0; % if zero, takes last. Seed is fixed  for random sampling
        targetcell = equalizegonogotrials(targetcell,0,0,randeq);% frnormalize
    end
    
    animalmodel{animalnum}.eqtrials = lmodel.eqtrials;
    animalmodel{animalnum}.randeq = randeq;
    animalmodel{animalnum}.type = lmodel.type;
    animalmodel{animalnum}.Xth = lmodel.Xth;
    animalmodel{animalnum}.hyperparams.gammarange = lmodel.hyperparams.gammarange;
    animalmodel{animalnum}.hyperparams.deltarange = lmodel.hyperparams.deltarange;
    animalmodel{animalnum}.hyperparams.cvn = lmodel.cvn;
    
    for gammai = 1:length(lmodel.hyperparams.gammarange)
        for deltai = 1:length(lmodel.hyperparams.deltarange)
            animalnum
            gammai
            lmodel.hyperparams.gamma = lmodel.hyperparams.gammarange(gammai); %0.2 0.6-0.9 %80ms: 0.8  per animal:0.2 - 80ms 0.8
            lmodel.hyperparams.delta = lmodel.hyperparams.deltarange(deltai); %0.2 %80ms: 0.2 per animal:0.2 - 80 ms 0.1
            
            %%%
            % gonig through all lags, getting vectors, perf, etc
            for i=1:nlags
                X = [cell2mat(cellfun(@(x) x.laAbs.go{i}',targetcell,'UniformOutput',0)')';...
                    cell2mat(cellfun(@(x) x.laAls.go{i}',targetcell,'UniformOutput',0)')'];
                Y = [zeros(size(targetcell{i}.laAbs.go{i},1),1);ones(size(targetcell{i}.laAls.go{i},1),1)];
                % alltrials*ncells*time
                alltraces = make_alltraces(X,targetcell,i,'go',lmodel);
                % do not manipulate X,Y befopre this point (proj depends on it)
                [vec_go{i},cnst_go{i},vec_n_go{i},~,cv_go{i},~,~,~,~] = fit_linear_classifier(X,Y,lmodel,alltraces,proj_vec_go);
                
                X = [cell2mat(cellfun(@(x) x.laAbs.nogo{i}',targetcell,'UniformOutput',0)')';...
                    cell2mat(cellfun(@(x) x.laAls.nogo{i}',targetcell,'UniformOutput',0)')'];
                Y = [zeros(size(targetcell{i}.laAbs.nogo{i},1),1);ones(size(targetcell{i}.laAls.nogo{i},1),1)];
                alltraces = make_alltraces(X,targetcell,i,'nogo',lmodel);
                [vec_nogo{i},cnst_nogo{i},vec_n_nogo{i},~,cv_nogo{i},~,~,~,~] = fit_linear_classifier(X,Y,lmodel,alltraces,proj_vec_nogo);
                
            end
            
            animalmodel{animalnum}.MeanTestEr(gammai,deltai) = mean([cellfun(@(x) x.test_er,cv_go),cellfun(@(x) x.test_er,cv_nogo)]);
            animalmodel{animalnum}.MeanTrainEr(gammai,deltai) = mean([cellfun(@(x) x.train_er,cv_go),cellfun(@(x) x.train_er,cv_nogo)]);
        end
    end
    
end
%%% save to animalfile
if savetofile
    cd(savepath);
    savename = sprintf('hyper_%slmodel25Th65ms_C%d_R%d_CV%d.mat',lmodel.exptype,lmodel.clean,lmodel.remove,lmodel.cvn);
    save(savename, 'animalmodel')
   
end
%% inspection: best hyper params per animal
% Here the minimum value of error is chosen, 


savepath = '/mnt/data/Mitra/figs/P2_L/bothdircombined/2020_02_24/LDA/';
savetofile = 1;

cd(savepath)
savename = 'hyper_FFlmodel25Th65ms_C1_R1_CV10.mat';

load(savename);



for animalnum=1:numanimals
    % plot test error as a function of gamma and delta
    figure;
    imagesc(animalmodel{animalnum}.hyperparams.gammarange,animalmodel{animalnum}.hyperparams.deltarange,...
        animalmodel{animalnum}.MeanTestEr)
    % find minimum test error
    animalmodel{animalnum}.mintester = min(min(animalmodel{animalnum}.MeanTestEr));
    animalmodel{animalnum}.mintesterdelta0 = min(animalmodel{animalnum}.MeanTestEr(:,1));
    
     [bestgammai,bestdeltai] = find(animalmodel{animalnum}.MeanTestEr == animalmodel{animalnum}.mintester,1);
     [bestgammaidelta0] = find(animalmodel{animalnum}.MeanTestEr(:,1) == animalmodel{animalnum}.mintesterdelta0,1);
    % best gamma and delta to produce minimum error
    animalmodel{animalnum}.bestgamma = animalmodel{animalnum}.hyperparams.gammarange(bestgammai);
    animalmodel{animalnum}.bestdelta = animalmodel{animalnum}.hyperparams.deltarange(bestdeltai);
    animalmodel{animalnum}.bestgammadelta0 = animalmodel{animalnum}.hyperparams.gammarange(bestgammaidelta0);
    % show best gamma, delta
    [animalnum,...
        animalmodel{animalnum}.bestgamma,...
        animalmodel{animalnum}.bestdelta,...
        animalmodel{animalnum}.bestgammadelta0]
    
end

if savetofile
    cd(savepath);
    save(savename, 'animalmodel','-append');       
end
