

%% params
% histogram plots
hist_yl = [0 0.04]; % ylimit
hist_edges = -100:10:200; % covered time in ms
hist_edges_bts = -100:2:200; % used for bootstraped distribution: finer bins - normalizaion should be pdf, when 2 bin types are used
mean_standarderr_coeff =1.96;% 1.96
nrep = 100;% 100
removenans = 1;
%%
% set V1 or LM cell
targetcell = eval(targetcellname);
% subselect cells for plots only
if subselectcells
    subselect_ind = eval(subselect_ind_phrase);
    fprintf('subselecting %d cells\n',numel(subselect_ind))
    targetcell = targetcell(subselect_ind);
    if showspikeshape
        allspikes = cell2mat(cellfun(@(x) x.spiketemplate,targetcell,'UniformOutput',0));
        figure;plot(allspikes(:,:),'k')
    end
end
%% preparing data

allgo = cell2mat(cellfun(@(x) x.smb.go, targetcell,'UniformOutput',0));
allnogo = cell2mat(cellfun(@(x) x.smb.nogo, targetcell,'UniformOutput',0));
all = cell2mat([cellfun(@(x) x.smb.nogo, targetcell,'UniformOutput',0),cellfun(@(x) x.smb.go, targetcell,'UniformOutput',0)]);

% standard error of the mean
standarderr_prct_go = nan(length(targetcell),length(targetcell{1}.laAbs.go));
standarderr_prct_nogo = nan(length(targetcell),length(targetcell{1}.laAbs.nogo));
for i=1:length(targetcell)
    % go
    for l=1:length(targetcell{i}.laAbs.go)
        mean_standarderr = nanstd(targetcell{i}.laAbs.go{l})/sqrt(numel(targetcell{i}.laAbs.go{l}));
        standarderr_prct_go(i,l) = 100 * mean_standarderr_coeff * mean_standarderr/nanmean(targetcell{i}.laAbs.go{l});
        if isinf(standarderr_prct_go(i,l))
            standarderr_prct_go(i,l) = nan;
        end
    end
    % nogo
    for l=1:length(targetcell{i}.laAbs.nogo)
        mean_standarderr = nanstd(targetcell{i}.laAbs.nogo{l})/sqrt(numel(targetcell{i}.laAbs.nogo{l}));
        standarderr_prct_nogo(i,l) = 100 * mean_standarderr_coeff * mean_standarderr/nanmean(targetcell{i}.laAbs.nogo{l});
        if isinf(standarderr_prct_nogo(i,l))
            standarderr_prct_nogo(i,l) = nan;
        end
    end
end


allmeanconfs = [reshape(standarderr_prct_nogo',1,[]),reshape(standarderr_prct_go',1,[])];
gomeanconfs = reshape(standarderr_prct_go',1,[]);
nogomeanconfs = reshape(standarderr_prct_nogo',1,[]);

% report the number of 'cells' in go/nogo that were ignificantly affected
% at least in 1 lag


all_s = all; allgo_s = allgo; allnogo_s = allnogo;
all_s(abs(all_s)<allmeanconfs) = nan;
allgo_s(abs(allgo_s)<gomeanconfs) = nan;
allnogo_s(abs(allnogo_s)<nogomeanconfs) = nan;

fprintf('%f%s of all go effects in %s are significant\n',...
    100*numel(find(~isnan(allgo_s)))/numel(find(~isnan(allgo))),'%',targetcellname(1:2))
fprintf('%f%s of all nogo effects in %s are significant\n',...
    100*numel(find(~isnan(allnogo_s)))/numel(find(~isnan(allnogo))),'%',targetcellname(1:2))
fprintf('%f%s of all effects in %s (go and nogo combined) are significant\n',...
    100*numel(find(~isnan(all_s)))/numel(find(~isnan(all))),'%',targetcellname(1:2))

% pie chart for sig positive and negative effects
sigpos = 100*numel(intersect(find(~isnan(all_s)),find(all_s>=0)))/numel(find(~isnan(all)));
signeg = 100*numel(intersect(find(~isnan(all_s)),find(all_s<0)))/numel(find(~isnan(all)));
nosigpos = 100*numel(intersect(find(isnan(all_s)),find(all>=0)))/numel(find(~isnan(all)));
nosigneg = 100*numel(intersect(find(isnan(all_s)),find(all<0)))/numel(find(~isnan(all)));
nonsig = 100*numel(find(isnan(all_s)))/numel(find(~isnan(all)));
figure;pie([sigpos,signeg,nonsig],[0 0 0],{num2str(sigpos),num2str(signeg),num2str(nonsig)});
legend({'pos','neg','nonsig'})
figure;pie([sigpos,nosigpos,signeg,nosigneg],[0 0 0 0],...
    {num2str(sigpos),num2str(nosigpos),num2str(signeg),num2str(nosigneg)});
legend({'pos','nosigpos','neg','nosigneg'})


% bootstraped differences: shuffled distribution: what would be the distribution if the silencing was not real
btstrp_prct_go = nan(length(targetcell),length(targetcell{1}.laAbs.go),nrep);
btstrp_prct_nogo = nan(length(targetcell),length(targetcell{1}.laAbs.nogo),nrep);
btstrp_standarderr_prct_go = nan(length(targetcell),length(targetcell{1}.laAbs.go),nrep);
btstrp_standarderr_prct_nogo = nan(length(targetcell),length(targetcell{1}.laAbs.nogo),nrep);

for i=1:length(targetcell)
    % go
    for l=1:length(targetcell{i}.laAbs.go)
        % for matching clean cell indices
        if ~isnan(targetcell{i}.smb.go(l))
            for r=1:nrep
                btstrpbs = targetcell{i}.laAbs.go{l};%datasample(targetcell{i}.laAbs.go{l},length(targetcell{i}.laAbs.go{l}));
                btstrpls = datasample(targetcell{i}.laAbs.go{l},length(targetcell{i}.laAls.go{l}));
                btstrp_prct_go(i,l,r) = 100*(nanmean(btstrpls) - nanmean(btstrpbs))/nanmean(btstrpbs);
                
                mean_standarderr = nanstd(btstrpbs)/sqrt(numel(btstrpbs));
                btstrp_standarderr_prct_go(i,l,r) = 100 * mean_standarderr_coeff * mean_standarderr/nanmean(btstrpbs);
                if isinf(btstrp_standarderr_prct_go(i,l,r))
                    btstrp_standarderr_prct_go(i,l,r) = nan;
                end
                
            end
        end
    end
    % nogo
    for l=1:length(targetcell{i}.laAbs.nogo)
        % for matching clean cell indices
        if ~isnan(targetcell{i}.smb.nogo(l))
            for r=1:nrep
                btstrpbs = targetcell{i}.laAbs.nogo{l};%datasample(targetcell{i}.laAbs.nogo{l},length(targetcell{i}.laAbs.nogo{l}));
                btstrpls = datasample(targetcell{i}.laAbs.nogo{l},length(targetcell{i}.laAls.nogo{l}));
                btstrp_prct_nogo(i,l,r) = 100*(nanmean(btstrpls) - nanmean(btstrpbs))/nanmean(btstrpbs);
               
                mean_standarderr = nanstd(btstrpbs)/sqrt(numel(btstrpbs));
                btstrp_standarderr_prct_nogo(i,l,r) = 100 * mean_standarderr_coeff * mean_standarderr/nanmean(btstrpbs);
                if isinf(btstrp_standarderr_prct_nogo(i,l,r))
                    btstrp_standarderr_prct_nogo(i,l,r) = nan;
                end
        
            end
        end
    end
end
% order is arbitrary
allbtstrp = [reshape(btstrp_prct_go,1,[]),reshape(btstrp_prct_nogo,1,[])];
%allbtstrp = cell2mat([cellfun(@(x) x.smb.nogo, targetcell,'UniformOutput',0),cellfun(@(x) x.smb.go, targetcell,'UniformOutput',0)]);

btstrp_allmeanconfs = [reshape(btstrp_standarderr_prct_go,1,[]),reshape(btstrp_standarderr_prct_nogo,1,[])]; 

allbtstrp_s = allbtstrp; 
allbtstrp_s(abs(allbtstrp_s)<btstrp_allmeanconfs) = nan;

fprintf('%f%s of all effects in %s bootstrapped (go and nogo combined) are significant\n',...
    100*numel(find(~isnan(allbtstrp_s)))/numel(find(~isnan(allbtstrp))),'%',targetcellname(1:2))

%%%%%%%% removing nans: only to be done after assigning significance: index
%%%%%%%% matching
if removenans % order is important
    all_s(isnan(all)) = [];
    allgo_s(isnan(allgo)) = [];
    allnogo_s(isnan(allnogo)) = [];
    allbtstrp(isnan(allbtstrp)) = [];
    allbtstrp_s(isnan(allbtstrp)) = [];
    all(isnan(all))=[];
    allgo(isnan(allgo))=[];
    allnogo(isnan(allnogo))=[];
    
end

%% making figures
figure;
%% s1:overlaid go/nogo
s1=subplot(3,2,1);
s1.YLim = hist_yl;
% go
hold on;
histogram(allgo,hist_edges,...
    'Normalization',globalplot.histnorm,'FaceColor','g','FaceAlpha',0.1,'EdgeAlpha',0);
hold on;
histogram(allgo,hist_edges,...
    'Normalization',globalplot.histnorm,'EdgeColor',[0 0.8 0],'EdgeAlpha',1,'DisplayStyle','stairs');
hold on;
line([nanmedian(allgo),...
    nanmedian(allgo)],...
    hist_yl,'Color',[0 0.8 0],'Linewidth',1.5);
% indicate med value
text(10,.9*hist_yl(2),sprintf('go median: %d%s',round(nanmedian(allgo)),'%'),'FontSize',7);


% nogo
hold on;
histogram(allnogo,hist_edges,...
    'Normalization',globalplot.histnorm,'FaceColor','r','FaceAlpha',0.1,'EdgeAlpha',0);
hold on;
histogram(allnogo,hist_edges,...
    'Normalization',globalplot.histnorm,'EdgeColor',[0.8 0 0],'EdgeAlpha',1,'DisplayStyle','stairs');
hold on;
line([nanmedian(allnogo),...
    nanmedian(allnogo)],...
    hist_yl,'Color',[0.8 0 0],'Linewidth',1.5);
% indicate med value
text(10,.8*hist_yl(2),sprintf('nogo median: %d%s',round(nanmedian(allnogo)),'%'),'FontSize',7);
% setting axis properties
box off
xlim([hist_edges(1) hist_edges(end)]);

% go/nogo difference
gonopval = ranksum(allgo,allnogo);
text(10,.7*hist_yl(2),sprintf('ranksum pval: %f',gonopval),'FontSize',7)


%% s2: go
s2=subplot(3,2,3);
s2.YLim = hist_yl;
hold on;
histogram(allgo,hist_edges,...
    'Normalization',globalplot.histnorm,'FaceColor','g','FaceAlpha',0.1,'EdgeAlpha',0);
hold on;
histogram(allgo,hist_edges,...
    'Normalization',globalplot.histnorm,'EdgeColor',[0 0.8 0],'EdgeAlpha',1,'DisplayStyle','stairs');
hold on;
line([nanmedian(allgo),...
    nanmedian(allgo)],...
    hist_yl,'Color',[0 0.8 0],'Linewidth',1.5);
% sig
hold on;
histogram(allgo_s,...
    hist_edges,'FaceColor',[0 1 0],'FaceAlpha',.8,'EdgeAlpha',0,'Normalization',globalplot.histnorm);

% indicate med value
line([nanmedian(allgo_s), nanmedian(allgo_s)],...
    hist_yl,'Color',[0 0.8 0],'Linewidth',1.5);
% text(0,.7*hist_yl(2),sprintf('median: %d%s',round(nanmedian(allgo_s)),'%'));

% setting axis properties
box off
xlim([hist_edges(1) hist_edges(end)])

%% s3: nogo
s3=subplot(3,2,5);
s3.YLim = hist_yl;
hold on;
histogram(allnogo,hist_edges,...
    'Normalization',globalplot.histnorm,'FaceColor','r','FaceAlpha',0.1,'EdgeAlpha',0);
hold on;
histogram(allnogo,hist_edges,...
    'Normalization',globalplot.histnorm,'EdgeColor',[0.8 0 0],'EdgeAlpha',1,'DisplayStyle','stairs');
hold on;
line([nanmedian(allnogo),...
    nanmedian(allnogo)],...
    hist_yl,'Color',[0.8 0 0],'Linewidth',1.5);
% sig
hold on;
histogram(allnogo_s,...
    hist_edges,'FaceColor',[1 0 0],'FaceAlpha',.8,'EdgeAlpha',0,'Normalization',globalplot.histnorm);

% indicate med value
line([nanmedian(allnogo_s), nanmedian(allnogo_s)],...
    hist_yl,'Color',[0.8 0 0],'Linewidth',1.5);
%text(0,.7*hist_yl(2),sprintf('median: %d%s',round(nanmedian(allnogo_s)),'%'));

% setting axis properties
box off
xlim([hist_edges(1) hist_edges(end)])


%% s4: pooled go/nogo
s4=subplot(3,2,2);
s4.YLim = hist_yl;
hold on;
histogram(all,...
    hist_edges,'FaceColor',[0 0 0],'FaceAlpha',0.8,'EdgeAlpha',0,'Normalization',globalplot.histnorm);
line([nanmedian(all),...
    nanmedian(all)],...
    hist_yl,'Color','k','Linewidth',1.5);
% indicate med value
%text(0,.8*hist_yl(2),sprintf('median: %d%s',round(nanmedian(all)),'%'));
% sig
hold on;
histogram(all_s,...
    hist_edges,'FaceColor',[0 0 0],'FaceAlpha',.9,'EdgeAlpha',0,'Normalization',globalplot.histnorm);

% indicate med value
line([nanmedian(all_s),nanmedian(all_s)],...
    hist_yl,'Color','k','Linewidth',1.5);

% indicate med value
text(10,.8*hist_yl(2),sprintf('median: %d%s',round(nanmedian(all)),'%'),'FontSize',7);

% setting axis properties
box off
xlim([hist_edges(1) hist_edges(end)])

%% s5: pooled go/nogo
s5=subplot(3,2,4);
s5.YLim = hist_yl;
hold on;
histogram(all,...
    hist_edges,'FaceColor',[0 0 1],'FaceAlpha',0.5,'EdgeAlpha',0,'Normalization',globalplot.histnorm);
line([nanmedian(all),...
    nanmedian(all)],...
    hist_yl,'Color','k','Linewidth',1.5);

% sig
hold on;
histogram(all_s,...
    hist_edges,'FaceColor',[0 0 1],'FaceAlpha',.9,'EdgeAlpha',0,'Normalization',globalplot.histnorm);

% add btstrp dist
hold on;
histogram(allbtstrp,...
    hist_edges_bts,'Normalization',globalplot.histnorm,'EdgeColor',[0 0 0],'EdgeAlpha',0.4,...
    'DisplayStyle','stairs','LineStyle','-');


% setting axis properties
box off
xlim([hist_edges(1) hist_edges(end)])


%% s6: bootstrap distribution for pooled go/nogo
s6=subplot(3,2,6);
s6.YLim = hist_yl;
hold on;
histogram(allbtstrp,...
    hist_edges_bts,'FaceColor',[0 0 0],'FaceAlpha',0.5,'EdgeAlpha',0,'Normalization',globalplot.histnorm);
line([nanmedian(allbtstrp),...
    nanmedian(allbtstrp)],...
    hist_yl,'Color','k','Linewidth',1.5);

% sig
hold on;
histogram(allbtstrp_s,...
    hist_edges_bts,'FaceColor',[0 0 0],'FaceAlpha',.9,'EdgeAlpha',0,'Normalization',globalplot.histnorm);




% setting axis properties
box off
xlim([hist_edges(1) hist_edges(end)])


%%

fprintf('sig go/nogo effect p value:%f\n',ranksum(allnogo_s,allgo_s))


