%% temporal change in pc direction, as the cosine similarity of pcs over time.
% + total variance

% [Latest version: completely from the lda file, no cell mat loaded.
% So clean and remove criteria is based on X_lda_new.m criteria]

warning('off')
exptype = 'FB';
binsizems = '150'; % default 150
xstyle = '1'; % default 1
doabs = 1; 
eqbsls = 4; 
numtrials = 10; % option to put an additional upper bound on the number of trials 
                % (doesn't change things much). nan or a number
nrep = 100;
totalvarnumcomp = nan; % nan to get total variance across all pcs, number to get first x pcs (plot both)

cd('/mnt/data/Mitra/figs/P2_L/bothdircombined/2020_02_24/LDA/version2')
load([exptype,'lmodel_25Th_',binsizems,'ms_X_style',xstyle,'_nrep',num2str(nrep),'.mat'])
includedanimals = [1 2 3 4 5 6]; % this is all feedback animals

%%%%% fit properties: default: trust-region
fo = fitoptions('Method','NonlinearLeastSquares','Lower',[-inf,-inf,0],'Upper',[inf,inf,inf],'StartPoint',[0.5,1,1000],'MaxIter',10000);
ftype = fittype('A*(exp(-x/tau)+B)','independent','x','options',fo);


timescalems = 65;
nlags = 8;
numcomp = 3;

nanimals = numel(includedanimals);
animal_min_trials = nan(nanimals,nrep);

% nlags*nanimals
totalbsvar.go = nan(8,nanimals);
totalbsvar.nogo = nan(8,nanimals);
totallsvar.go = nan(8,nanimals);
totallsvar.nogo = nan(8,nanimals);

lsexpvar.go = nan(8,nanimals);
lsexpvar.nogo = nan(8,nanimals);
bsexpvar.go = nan(8,nanimals);
bsexpvar.nogo = nan(8,nanimals);


bsvarspec.go = cell(8,nanimals);
bsvarspec.nogo = cell(8,nanimals);
lsvarspec.go = cell(8,nanimals);
lsvarspec.nogo = cell(8,nanimals);


for lsbs = 1:2
    
    gopcang = struct; gopcang.rep=cell(1,nrep);
    nogopcang = struct; nogopcang.rep=cell(1,nrep);
    
    for rep = 1:nrep
        
        gopcang.rep{rep}.t = cell(1,nlags);
        nogopcang.rep{rep}.t = cell(1,nlags);
        for t=1:8
            gopcang.rep{rep}.t{t} = nan(1,numcomp);
            nogopcang.rep{rep}.t{t} = nan(1,numcomp);
        end
        
        
        for animalnum = includedanimals
            
            if eqbsls == 1
                mintrials = min([cellfun(@(x) numel(find(x==0)),animalmodel{animalnum}.lmodel.rep{rep}.part{1}.go.Y),...
                    cellfun(@(x) numel(find(x==0)),animalmodel{animalnum}.lmodel.rep{rep}.part{2}.go.Y),...
                    cellfun(@(x) numel(find(x==0)),animalmodel{animalnum}.lmodel.rep{rep}.part{1}.nogo.Y),...
                    cellfun(@(x) numel(find(x==0)),animalmodel{animalnum}.lmodel.rep{rep}.part{2}.nogo.Y),...
                    cellfun(@(x) numel(find(x==1)),animalmodel{animalnum}.lmodel.rep{rep}.part{1}.go.Y),...
                    cellfun(@(x) numel(find(x==1)),animalmodel{animalnum}.lmodel.rep{rep}.part{2}.go.Y),...
                    cellfun(@(x) numel(find(x==1)),animalmodel{animalnum}.lmodel.rep{rep}.part{1}.nogo.Y),...
                    cellfun(@(x) numel(find(x==1)),animalmodel{animalnum}.lmodel.rep{rep}.part{2}.nogo.Y)]);
            elseif eqbsls == 2 || eqbsls == 3
                mintrials_bs = min([cellfun(@(x) numel(find(x==0)),animalmodel{animalnum}.lmodel.rep{rep}.part{1}.go.Y),...
                    cellfun(@(x) numel(find(x==0)),animalmodel{animalnum}.lmodel.rep{rep}.part{2}.go.Y),...
                    cellfun(@(x) numel(find(x==0)),animalmodel{animalnum}.lmodel.rep{rep}.part{1}.nogo.Y),...
                    cellfun(@(x) numel(find(x==0)),animalmodel{animalnum}.lmodel.rep{rep}.part{2}.nogo.Y)]);
                mintrials_ls = min([cellfun(@(x) numel(find(x==1)),animalmodel{animalnum}.lmodel.rep{rep}.part{1}.go.Y),...
                    cellfun(@(x) numel(find(x==1)),animalmodel{animalnum}.lmodel.rep{rep}.part{2}.go.Y),...
                    cellfun(@(x) numel(find(x==1)),animalmodel{animalnum}.lmodel.rep{rep}.part{1}.nogo.Y),...
                    cellfun(@(x) numel(find(x==1)),animalmodel{animalnum}.lmodel.rep{rep}.part{2}.nogo.Y)]);
                
                if lsbs == 1
                    mintrials = mintrials_ls;
                elseif lsbs == 2
                    mintrials = mintrials_bs;
                end
                
                if eqbsls == 3
                    mintrials = min(numtrials,mintrials);
                end
            elseif eqbsls == 4
                mintrials = nan(2,8);
                for lag = 1:8
                    
                        mintrials(1,lag) = min([numel(find(animalmodel{animalnum}.lmodel.rep{rep}.part{1}.go.Y{lag}==0)),...                      
                        numel(find(animalmodel{animalnum}.lmodel.rep{rep}.part{1}.nogo.Y{lag}==0)),...                 
                        numel(find(animalmodel{animalnum}.lmodel.rep{rep}.part{1}.go.Y{lag}==1)),...                      
                        numel(find(animalmodel{animalnum}.lmodel.rep{rep}.part{1}.nogo.Y{lag}==1))]);
                        
                        mintrials(2,lag) = min([numel(find(animalmodel{animalnum}.lmodel.rep{rep}.part{2}.go.Y{lag}==0)),...                      
                        numel(find(animalmodel{animalnum}.lmodel.rep{rep}.part{2}.nogo.Y{lag}==0)),...                 
                        numel(find(animalmodel{animalnum}.lmodel.rep{rep}.part{2}.go.Y{lag}==1)),...                      
                        numel(find(animalmodel{animalnum}.lmodel.rep{rep}.part{2}.nogo.Y{lag}==1))]);
                   
                end
            end
            
            
            for t = 0:7
                for  lag =  1:nlags-t
                    
                    population_go = struct; population_go.part = cell(1,2);
                    population_next_go = struct; population_next_go.part = cell(1,2);
                    population_nogo = struct; population_nogo.part = cell(1,2);
                    population_next_nogo = struct; population_next_nogo.part = cell(1,2);
                    
                    pcs_go = struct;
                    pcs_next_go = struct;
                    pcs_nogo = struct;
                    pcs_next_nogo = struct;
                    
                    pcs_go.part = cell(1,2);
                    pcs_next_go.part = cell(1,2);
                    pcs_nogo.part = cell(1,2);
                    pcs_next_nogo.part = cell(1,2);
                    
                    
                    for p = 1:2
   
                        
                        if lsbs == 2 % bs
                            
                            population_go.part{p} = animalmodel{animalnum}.lmodel.rep{rep}.part{p}.go.X{lag}...
                                (find(animalmodel{animalnum}.lmodel.rep{rep}.part{p}.go.Y{lag}==0),:);
                            population_next_go.part{p} = animalmodel{animalnum}.lmodel.rep{rep}.part{p}.go.X{lag+t}...
                                (find(animalmodel{animalnum}.lmodel.rep{rep}.part{p}.go.Y{lag+t}==0),:);
                            
                            population_nogo.part{p} = animalmodel{animalnum}.lmodel.rep{rep}.part{p}.nogo.X{lag}...
                                (find(animalmodel{animalnum}.lmodel.rep{rep}.part{p}.nogo.Y{lag}==0),:);
                            population_next_nogo.part{p} = animalmodel{animalnum}.lmodel.rep{rep}.part{p}.nogo.X{lag+t}...
                                (find(animalmodel{animalnum}.lmodel.rep{rep}.part{p}.nogo.Y{lag+t}==0),:);
                            
                        elseif lsbs == 1% ls
                            
                            population_go.part{p} = animalmodel{animalnum}.lmodel.rep{rep}.part{p}.go.X{lag}...
                                (find(animalmodel{animalnum}.lmodel.rep{rep}.part{p}.go.Y{lag}==1),:);
                            population_next_go.part{p} = animalmodel{animalnum}.lmodel.rep{rep}.part{p}.go.X{lag+t}...
                                (find(animalmodel{animalnum}.lmodel.rep{rep}.part{p}.go.Y{lag+t}==1),:);
                            
                            population_nogo.part{p} = animalmodel{animalnum}.lmodel.rep{rep}.part{p}.nogo.X{lag}...
                                (find(animalmodel{animalnum}.lmodel.rep{rep}.part{p}.nogo.Y{lag}==1),:);
                            population_next_nogo.part{p} = animalmodel{animalnum}.lmodel.rep{rep}.part{p}.nogo.X{lag+t}...
                                (find(animalmodel{animalnum}.lmodel.rep{rep}.part{p}.nogo.Y{lag+t}==1),:);
                            
                        end
                        
                        
                        
                        % if cells have been cleaned, remove nan before pca: (if indss not used)                        
                        population_go.part{p}(:,find(isnan(nanmean(population_go.part{p},1)))) = 0;
                        population_next_go.part{p}(:,find(isnan(nanmean(population_next_go.part{p},1)))) = 0;
                        population_nogo.part{p}(:,find(isnan(nanmean(population_nogo.part{p},1)))) = 0;
                        population_next_nogo.part{p}(:,find(isnan(nanmean(population_next_nogo.part{p},1)))) = 0;
                        
                        
                        %quick per animal and per lag equalizing of bs and ls number of trials
                        if eqbsls > 0 && eqbsls < 4
                            
                            population_go.part{p}((mintrials+1):end,:) = [];
                            population_next_go.part{p}((mintrials+1):end,:) = [];
                            population_nogo.part{p}((mintrials+1):end,:) = [];
                            population_next_nogo.part{p}((mintrials+1):end,:) = [];
                        elseif eqbsls == 4
                            population_go.part{p}((mintrials(p,lag)+1):end,:) = [];
                            population_next_go.part{p}((mintrials(p,lag+t)+1):end,:) = [];
                            population_nogo.part{p}((mintrials(p,lag)+1):end,:) = [];
                            population_next_nogo.part{p}((mintrials(p,lag+t)+1):end,:) = [];
                           
                        end
                        
                        if ~isnan(numtrials)                        
                             population_go.part{p}((numtrials+1):end,:) = [];
                             population_next_go.part{p}((numtrials+1):end,:) = [];
                             population_nogo.part{p}((numtrials+1):end,:) = [];
                             population_next_nogo.part{p}((numtrials+1):end,:) = [];
                        end

                        
                        [pcs_go.part{p},~,~,~,~]= pca(population_go.part{p},'Centered',1);
                        %  varsbsgo(end+1,1:length(var)) = var;
                        [pcs_next_go.part{p},~,~,~,~]= pca(population_next_go.part{p},'Centered',1);
                        %  varslsgo(end+1,1:length(var)) = var;
                        [pcs_nogo.part{p},~,~,~,~]= pca(population_nogo.part{p},'Centered',1);
                        % varsbsnogo(end+1,1:length(var)) = var;
                        [pcs_next_nogo.part{p},~,~,~,~]= pca(population_next_nogo.part{p},'Centered',1);
                        % varslsnogo(end+1,1:length(var)) = var;
                        
                        
                    end % part
                    
                    gopcang.rep{rep}.t{t+1}(end+1,:) = nan(1,numcomp);
                    nogopcang.rep{rep}.t{t+1}(end+1,:) = nan(1,numcomp);
                    
                    % checking the processed time windows have at least numcomp
                    % pcs. 
                    
                     if (size(pcs_go.part{1},2)>=numcomp && size(pcs_go.part{2},2)>=numcomp) && ...
                             (size(pcs_next_go.part{1},2)>=numcomp && size(pcs_next_go.part{2},2)>=numcomp) && ...
                             (size(pcs_nogo.part{1},2)>=numcomp && size(pcs_nogo.part{2},2)>=numcomp) && ...
                             (size(pcs_next_nogo.part{1},2)>=numcomp && size(pcs_next_nogo.part{2},2)>=numcomp)
                         
                    
                                              
                        % calculate variance at each laser lag, _next
                        % variables are not relevant. so t = 0
                        % also both parts are combined, so rep does not
                        % matter
                        
                        if t == 0 && rep == 1
                           
                            if lsbs == 1
                                [totallsvar.go(lag,animalnum),totallsvar.nogo(lag,animalnum)] = calculate_total_var(population_go,population_nogo,totalvarnumcomp);
                                [lsvarspec.go{lag,animalnum},lsvarspec.nogo{lag,animalnum}] = calculate_spec_var(population_go,population_nogo);
                                [lsexpvar.go(lag,animalnum),lsexpvar.nogo(lag,animalnum)] = calculate_var_exp(population_go,population_nogo,numcomp);
                            elseif lsbs == 2
                                [totalbsvar.go(lag,animalnum),totalbsvar.nogo(lag,animalnum)] = calculate_total_var(population_go,population_nogo,totalvarnumcomp);
                                [bsvarspec.go{lag,animalnum},bsvarspec.nogo{lag,animalnum}] = calculate_spec_var(population_go,population_nogo);
                                [bsexpvar.go(lag,animalnum),bsexpvar.nogo(lag,animalnum)] = calculate_var_exp(population_go,population_nogo,numcomp);
                            end
                           
                        end
                   
                        for i=1:numcomp
                            try
                                half1 = (pcs_go.part{1}(:,i)'*pcs_next_go.part{2}(:,i));
                            catch
                                half1 = nan;
                            end
                            try
                                half2 = (pcs_go.part{2}(:,i)'*pcs_next_go.part{1}(:,i));
                            catch
                                half2 = nan;
                            end
                      
                            gopcang.rep{rep}.t{t+1}(end,i) = nanmean([half1,half2]);
                           
                        end
                        for  i=1:numcomp
                            try
                                half1 = (pcs_nogo.part{1}(:,i)'*pcs_next_nogo.part{2}(:,i));
                            catch
                                half1 = nan;
                            end
                            try 
                                half2 = (pcs_nogo.part{2}(:,i)'*pcs_next_nogo.part{1}(:,i));
                            catch
                                half2 = nan;
                            end
                            
                            nogopcang.rep{rep}.t{t+1}(end,i) =  nanmean([half1,half2]);
                           
                            
                        end
                    end
            
                    
                end % lag
                
            end % t
            
            
        end % animal
        
        
        if doabs
            for t=1:8
                gopcang.rep{rep}.t{t} = abs(gopcang.rep{rep}.t{t});
                nogopcang.rep{rep}.t{t} = abs(nogopcang.rep{rep}.t{t});
            end
        end
        
        
    end % rep
    
    %%%% average reps now
    gopcangAv = cell(1,nlags);
    nogopcangAv = cell(1,nlags);
    for lag = 1:nlags
        gopcangAv{lag} = zeros(size(gopcang.rep{1}.t{lag}));
        nogopcangAv{lag} = zeros(size(nogopcang.rep{1}.t{lag}));
        for rep = 1:nrep
            gopcangAv{lag} =gopcangAv{lag} + gopcang.rep{rep}.t{lag};
            nogopcangAv{lag} =nogopcangAv{lag} + nogopcang.rep{rep}.t{lag};
        end
        %%% abs could be added here if not before
        gopcangAv{lag} = abs(gopcangAv{lag}/nrep);
        nogopcangAv{lag} = abs(nogopcangAv{lag}/nrep);
    end
    
    % replace the structs with rep averages
    clear gopcang nogopcang
    gopcang = gopcangAv;
    nogopcang = nogopcangAv;
    
    %%%%%%%%% figures
    
    gopcang_i = gopcang;
    nogopcang_i = nogopcang;
    
    figure
    for pcn = 1:numcomp
        
        gopcang = gopcang_i;
        nogopcang = nogopcang_i;
        
        gomeans = nan(1,7);
        nogomeans = nan(1,7);
        for t=1:8
            s=subplot(1,numcomp,pcn);
            errorbar(timescalems*(t-1),nanmean(gopcang{t}(:,pcn)),2*nanstd(gopcang{t}(:,pcn))/sqrt(numel(find(~isnan(gopcang{t}(:,pcn))))),'g');
            %bc=bootci(100,@nanmean,gopcang{t}(:,pcn));
            %errorbar(timescalems*(t-1),nanmean(gopcang{t}(:,pcn)),abs(bc(1)-nanmean(gopcang{t}(:,pcn))),abs(bc(2)-nanmean(gopcang{t}(:,pcn))),'g');
            hold on;
            errorbar(timescalems*(t-1),nanmean(nogopcang{t}(:,pcn)),2*nanstd(nogopcang{t}(:,pcn))/sqrt(numel(find(~isnan(nogopcang{t}(:,pcn))))),'r');
            %bc=bootci(100,@nanmean,nogopcang{t}(:,pcn));
            %errorbar(timescalems*(t-1),nanmean(nogopcang{t}(:,pcn)),abs(bc(1)-nanmean(nogopcang{t}(:,pcn))),abs(bc(2)-nanmean(nogopcang{t}(:,pcn))),'r');
            hold on;
            scatter(timescalems*(t-1),nanmean(gopcang{t}(:,pcn)),'g.')
            hold on;
            scatter(timescalems*(t-1),nanmean(nogopcang{t}(:,pcn)),'r.')
            ylim([-0,0.8])
            %ylim([-0,30])
            s.Title.String = ['pc',num2str(pcn)];
            
            gomeans(t) = nanmean(gopcang{t}(:,pcn));
            nogomeans(t) = nanmean(nogopcang{t}(:,pcn));
            
        end
        
        
      
        %%%%%
        if pcn >= 1
            
            xgo = [zeros(size(gopcang{1}(:,pcn)));1+zeros(size(gopcang{2}(:,pcn)));2+zeros(size(gopcang{3}(:,pcn)));...
                3+zeros(size(gopcang{4}(:,pcn)));4+zeros(size(gopcang{5}(:,pcn)));5+zeros(size(gopcang{6}(:,pcn)));...
                6+zeros(size(gopcang{7}(:,pcn)));7+zeros(size(gopcang{8}(:,pcn)))];
            ygo = [gopcang{1}(:,pcn);gopcang{2}(:,pcn);gopcang{3}(:,pcn);gopcang{4}(:,pcn);...
                gopcang{5}(:,pcn);gopcang{6}(:,pcn);gopcang{7}(:,pcn);gopcang{8}(:,pcn)];
            xnogo = [zeros(size(nogopcang{1}(:,pcn)));1+zeros(size(nogopcang{2}(:,pcn)));2+zeros(size(nogopcang{3}(:,pcn)));...
                3+zeros(size(nogopcang{4}(:,pcn)));4+zeros(size(nogopcang{5}(:,pcn)));5+zeros(size(nogopcang{6}(:,pcn)));...
                6+zeros(size(nogopcang{7}(:,pcn)));7+zeros(size(nogopcang{8}(:,pcn)))];
            ynogo = [nogopcang{1}(:,pcn);nogopcang{2}(:,pcn);nogopcang{3}(:,pcn);nogopcang{4}(:,pcn);...
                nogopcang{5}(:,pcn);nogopcang{6}(:,pcn);nogopcang{7}(:,pcn);nogopcang{8}(:,pcn)];
            
            novals=isnan(ygo);
            xgo(novals)=[];ygo(novals)=[];
            
            novals=isnan(ynogo);
            xnogo(novals)=[];ynogo(novals)=[];
            
            xgo = timescalems*xgo;
            xnogo = timescalems*xnogo;
            
            try
                
                [mdlgo,gof] = fit(xgo,ygo,ftype);
                [mdlnogo,gof] = fit(xnogo,ynogo,ftype);
                
                subplot(1,numcomp,pcn);
                hold on;
                plot(unique(xgo),mdlgo.A*(exp(-(unique(xgo))/mdlgo.tau)+mdlgo.B),'-','Color','g')
                hold on;
                plot(unique(xnogo),mdlnogo.A*(exp(-(unique(xnogo))/mdlnogo.tau)+mdlnogo.B),'-','Color','r')
                ypgo = predint(mdlgo,unique(xgo),0.95,'functional','on');
                ypnogo = predint(mdlnogo,unique(xnogo),0.95,'functional','on');
                patch([unique(xgo)',fliplr(unique(xgo)')],...
                    [ypgo(:,1)',fliplr(ypgo(:,2)')],'g',...
                    'EdgeAlpha',0,'FaceAlpha',0.2);
                hold on;
                patch([unique(xnogo)',fliplr(unique(xnogo)')],...
                    [ypnogo(:,1)',fliplr(ypnogo(:,2)')],'r',...
                    'EdgeAlpha',0,'FaceAlpha',0.2);
                hold on;
                %pause(1)
            catch
                disp('no fit')
            end
        end
        xlim(timescalems*[-1 9])
        
        hold on;text(1, 0.05,sprintf('tau = %f, A = %f, B = %f',mdlgo.tau,mdlgo.A,mdlgo.B),'Color','g','FontSize',6);
        hold on;text(1, 0.1,sprintf('tau = %f, A = %f, B = %f',mdlnogo.tau,mdlnogo.A,mdlnogo.B),'Color','r','FontSize',6);
    end
    title(['lsbs = ',num2str(lsbs)])
    
    pause(1);
end


%% total variance plots
% totalbsvar.nogo type variables are nlags*nanimals.
% each entry is the sum of variance along all pcs
% the distributions are over both lags and animals
try
figure;
% go

bar(1,nanmean(reshape(totalbsvar.go,1,[])),'g','FaceAlpha',0.5,'EdgeAlpha',0);hold on;
errorbar(1,nanmean(reshape(totalbsvar.go,1,[])),...
    2*nanstd(reshape(totalbsvar.go,1,[]))/sqrt(numel(find(~isnan(totalbsvar.go)))),'g');hold on;

bar(2,nanmean(reshape(totallsvar.go,1,[])),'g','FaceAlpha',0.2,'EdgeAlpha',0);hold on;
errorbar(2,nanmean(reshape(totallsvar.go,1,[])),...
    2*nanstd(reshape(totallsvar.go,1,[]))/sqrt(numel(find(~isnan(totallsvar.go)))),'g');hold on;
text(1,45,['ranksum p  = ',num2str(ranksum(reshape(totalbsvar.go,1,[]),reshape(totallsvar.go,1,[])))])
text(1,47,['singrank p  = ',num2str(signrank(reshape(totalbsvar.go,1,[]),reshape(totallsvar.go,1,[])))])
%ylim([0 60]);
%xticks([1,2]);
%xticklabels({'bs-go','ls-go'})

% nogo
hold on;
bar(4,nanmean(reshape(totalbsvar.nogo,1,[])),'r','FaceAlpha',0.5,'EdgeAlpha',0);hold on;
errorbar(4,nanmean(reshape(totalbsvar.nogo,1,[])),...
    2*nanstd(reshape(totalbsvar.nogo,1,[]))/sqrt(numel(find(~isnan(totalbsvar.nogo)))),'r');hold on;

bar(5,nanmean(reshape(totallsvar.nogo,1,[])),'r','FaceAlpha',0.2,'EdgeAlpha',0);hold on;
errorbar(5,nanmean(reshape(totallsvar.nogo,1,[])),...
    2*nanstd(reshape(totallsvar.nogo,1,[]))/sqrt(numel(find(~isnan(totallsvar.nogo)))),'r');hold on;
text(4,45,['ranksum p  = ',num2str(ranksum(reshape(totalbsvar.nogo,1,[]),reshape(totallsvar.nogo,1,[])))])
text(4,47,['signrank p  = ',num2str(signrank(reshape(totalbsvar.nogo,1,[]),reshape(totallsvar.nogo,1,[])))])
ylim([0 60])

xticks([1,2,4,5]);
xticklabels({'bs-go','ls-go','bs-nogo','ls-nogo'});

text(3,55,['bs go nogo ranksum p  = ',num2str(ranksum(reshape(totalbsvar.go,1,[]),reshape(totalbsvar.nogo,1,[])))])
text(3,57,['bs go nogo signrank p  = ',num2str(signrank(reshape(totalbsvar.go,1,[]),reshape(totalbsvar.nogo,1,[])))])

text(3,53,['ls go nogo ranksum p  = ',num2str(ranksum(reshape(totallsvar.go,1,[]),reshape(totallsvar.nogo,1,[])))])
text(3,51,['ls go nogo signrank p  = ',num2str(signrank(reshape(totallsvar.go,1,[]),reshape(totallsvar.nogo,1,[])))])
%% explained variance: calculated in bs trials (by numcomp pcs compared to overal variance)
expm = nanmean([reshape(bsexpvar.go,1,[]),reshape(bsexpvar.nogo,1,[])]);
expstd = nanstd([reshape(bsexpvar.go,1,[]),reshape(bsexpvar.nogo,1,[])]);
sprintf('percentage explained variance by %d pcs = %f (%f) [mean(std)]',numcomp,expm,expstd)
%% pool across animals and bins: only bins with at least 3 pcs were processed, 
% and only the forst 3 pcs are used for eigen spectrum 
% (otherwise, pooling across animals with different number of dimensions not trivial)
numspec_pc = 3;
% remove animals/timewinsows that have less than numspec_pc pcs
bsvarspecgo_l=reshape(bsvarspec.go,1,[]);
lsvarspecgo_l=reshape(lsvarspec.go,1,[]);
bsvarspecnogo_l=reshape(bsvarspec.nogo,1,[]);
lsvarspecnogo_l=reshape(lsvarspec.nogo,1,[]);

bsvarspecgo_l(find(cellfun(@(x) numel(x), bsvarspecgo_l)==0)) = [];
lsvarspecgo_l(find(cellfun(@(x) numel(x), lsvarspecgo_l)==0)) = [];
bsvarspecnogo_l(find(cellfun(@(x) numel(x), bsvarspecnogo_l)==0)) = [];
lsvarspecnogo_l(find(cellfun(@(x) numel(x), lsvarspecnogo_l)==0)) = [];

% check here
bsgo = cell2mat(cellfun(@(x) x(1:numspec_pc)/nansum(x),bsvarspecgo_l,'UniformOutput',0));
lsgo = cell2mat(cellfun(@(x) x(1:numspec_pc)/nansum(x),lsvarspecgo_l,'UniformOutput',0));
bsnogo = cell2mat(cellfun(@(x) x(1:numspec_pc)/nansum(x),bsvarspecnogo_l,'UniformOutput',0));
lsnogo = cell2mat(cellfun(@(x) x(1:numspec_pc)/nansum(x),lsvarspecnogo_l,'UniformOutput',0));

figure;subplot(1,2,1);
plwdt= .05;
plot([1:numspec_pc]+plwdt,nanmean(bsgo,2),'g')
hold on;plot([1:numspec_pc]-plwdt,nanmean(lsgo,2),'b')
hold on;errorbar([1:numspec_pc]+plwdt,nanmean(bsgo,2),2*nanstd(bsgo')/sqrt(size(bsgo,2)),'g')
hold on;errorbar([1:numspec_pc]-plwdt,nanmean(lsgo,2),2*nanstd(lsgo')/sqrt(size(lsgo,2)),'b')
xlabel('pc number');ylabel('variance')
xlim([-0.5,4.5])
ylim([0 1.2])
text(1,0.8,['signrank p  = ',num2str(signrank(bsgo(1,:),lsgo(1,:)))])
text(2,0.5,['signrank p  = ',num2str(signrank(bsgo(2,:),lsgo(2,:)))])
text(3,0.3,['signrank p  = ',num2str(signrank(bsgo(3,:),lsgo(3,:)))])

subplot(1,2,2);
plot([1:numspec_pc]+plwdt,nanmean(bsnogo,2),'r')
hold on;plot([1:numspec_pc]-plwdt,nanmean(lsnogo,2),'b')
hold on;errorbar([1:numspec_pc]+plwdt,nanmean(bsnogo,2),2*nanstd(bsnogo')/sqrt(size(bsnogo,2)),'r')
hold on;errorbar([1:numspec_pc]-plwdt,nanmean(lsnogo,2),2*nanstd(lsnogo')/sqrt(size(lsnogo,2)),'b')
xlabel('pc number');ylabel('variance')
xlim([-0.5,4.5])
ylim([0 1.2])
text(1,1,['signrank p  = ',num2str(signrank(bsnogo(1,:),lsnogo(1,:)))])
text(2,0.5,['signrank p  = ',num2str(signrank(bsnogo(2,:),lsnogo(2,:)))])
text(3,0.3,['signrank p  = ',num2str(signrank(bsnogo(3,:),lsnogo(3,:)))])
catch
end

%%
function [tvargo tvarnogo] = calculate_total_var(population_go,population_nogo,totalvarnumcomp)
% go
% numalltrials * numcells
allgo = [population_go.part{1};population_go.part{2}];
[~,~,vargo]= pca(allgo,'Centered',1);
% var go is one number per pc. Sum of all its elements is total var
if isnan(totalvarnumcomp)
    tvargo = sum(vargo);
else
    tvargo = sum(vargo(1:totalvarnumcomp));
end
    


%nogo
% numalltrials * numcells
allnogo = [population_nogo.part{1};population_nogo.part{2}];
[~,~,varnogo]= pca(allnogo,'Centered',1);
% var go is one number per pc. Sum of all its elements is total var
if isnan(totalvarnumcomp)
    tvarnogo = sum(varnogo);
else
    tvarnogo = sum(varnogo(1:totalvarnumcomp));
end

end

function [tvargo tvarnogo] = calculate_spec_var(population_go,population_nogo)
% go

% numalltrials * numcells
allgo = [population_go.part{1};population_go.part{2}];
[~,~,vargo]= pca(allgo,'Centered',1);
% var go is one number per pc
% only lags with at least 3 pcs are included.  (either due to numcells or num trials)
tvargo = vargo;

%nogo

% numalltrials * numcells
allnogo = [population_nogo.part{1};population_nogo.part{2}];
[~,~,varnogo]= pca(allnogo,'Centered',1);
tvarnogo = varnogo;

end

function [tvargo tvarnogo] = calculate_var_exp(population_go,population_nogo,totalvarnumcomp)
% go
% numalltrials * numcells
allgo = [population_go.part{1};population_go.part{2}];
[~,~,vargo]= pca(allgo,'Centered',1);
tvargo = sum(vargo(1:totalvarnumcomp))/sum(vargo);

%nogo
% numalltrials * numcells
allnogo = [population_nogo.part{1};population_nogo.part{2}];
[~,~,varnogo]= pca(allnogo,'Centered',1);
tvarnogo = sum(varnogo(1:totalvarnumcomp))/sum(varnogo);


end