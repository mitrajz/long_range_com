function Aggregate_BaselineGovsNogoPlots(V1cells,LMcells,params,plotstodo,ntrial_th)
%%
stimlengthinms = 500; %500
zoomedin=0;
doerbars = 1;
nrepbtstrp = 100; % def 1000
alphabtstrp = 5;
clm=0.5;
indv_sorting = 1;
sortingtype = 'go'; % only counts when indv_sorting = 0
% which plots to do
smoothingspan = 1; % only used for traces


timeoffset = .5*params.edgestep_pl/(params.samplingf/1000);


 selindex=1; 
%% raw traces og LM and V1 (all units) in go vs. nogo
% set V1 or LM cell
xpoints = ((-params.Window):1*params.edgestep_pl:(params.Window-1))/(params.samplingf/1000);
targetcellnamelist = {'V1cells','LMcells'};

for areaind = 1:2
    
    targetcellname = targetcellnamelist{areaind};
    targetcell = eval(targetcellname);
    
    zerotrialsinonecond = find(cellfun(@(x) min(size(x.nbs.go,1),size(x.nbs.nogo,1)),targetcell)<ntrial_th);
    targetcell(zerotrialsinonecond)=[];
        
    allgo=cell2mat(transpose(cellfun(@(x) smooth(x.nbsAv.go,smoothingspan)' ,targetcell,'UniformOutput',0)));
    allnogo=cell2mat(transpose(cellfun(@(x) smooth(x.nbsAv.nogo,smoothingspan)' ,targetcell,'UniformOutput',0)));
    
    %% traces
    if plotstodo.trace
        figure;plot(timeoffset + xpoints,1000*nanmean(allgo,1)/(params.edgestep_pl/30),'Color',[0 0.8 0],'LineWidth',1.6);
        hold on;plot(timeoffset + xpoints,1000*nanmean(allnogo,1)/(params.edgestep_pl/30),'Color',[1 0 0],'LineWidth',1.6);
        
        if doerbars
            allgo_rand_mean=nan(nrepbtstrp,size(allgo,2));
            allnogo_rand_mean=nan(nrepbtstrp,size(allnogo,2));
            for randsamp = 1:nrepbtstrp
                allgo_rand = cell2mat(transpose(cellfun(@(x) smooth(nanmean(x.nbs.go(randi(size(x.nbs.go,1),1,size(x.nbs.go,1)),:),1) , smoothingspan)',...
                    targetcell,'UniformOutput',0)));
                allnogo_rand = cell2mat(transpose(cellfun(@(x) smooth(nanmean(x.nbs.nogo(randi(size(x.nbs.nogo,1),1,size(x.nbs.nogo,1)),:),1) , smoothingspan)',...
                    targetcell,'UniformOutput',0)));
                allgo_rand_mean(randsamp,:) = nanmean(allgo_rand,1)/(params.edgestep_pl/30);
                allnogo_rand_mean(randsamp,:) = nanmean(allnogo_rand,1)/(params.edgestep_pl/30);
            end
            hold on;patch([timeoffset + xpoints,fliplr(timeoffset + xpoints)],...
                [1000*prctile(allgo_rand_mean,100 - (alphabtstrp/2)),fliplr(1000*prctile(allgo_rand_mean,(alphabtstrp/2)))],[0 0.8 0],'FaceAlpha',0.2,'EdgeAlpha',0);
            hold on;patch([timeoffset + xpoints,fliplr(timeoffset + xpoints)],...
                [1000*prctile(allnogo_rand_mean,100 - (alphabtstrp/2)),fliplr(1000*prctile(allnogo_rand_mean,(alphabtstrp/2)))],[1 0 0],'FaceAlpha',0.2,'EdgeAlpha',0);
            
            
        end
        
        xlabel('Time (ms)');ylabel('Firing rate (Hz)');
        xlim([-100 650]);ylim([0 25]);
        if zoomedin
            ylim([0 25]);
        end
        axx = gca;
        hold on; patch([0 stimlengthinms stimlengthinms 0],[axx.YLim(1) axx.YLim(1) axx.YLim(2) axx.YLim(2)],[0.7 0.7 0.7],...
            'FaceAlpha',0.3,'EdgeColor','none')
        axx.XTick=0:100:650;
        axx.XTickLabel={'0','','200','','400','','600',''};
        if zoomedin
            xlim([0 200]);
            axx.XTick=0:50:200;
            axx.XTickLabel={'0','','100','','200'};
        end
        axx.YTick=0:5:20;
        axx.YTickLabel={'0','5','10','15','20','25'};
        set(gcf,'Color','w')
        legend({'go - correct','no go - correct'})
        legend('boxoff')
        box(axx,'off');
        if zoomedin
            legend('hide');
        end
    end
    %% imagesc: sorted off go or nogo responsivity
    if plotstodo.raster
        [~,sorti_go] = sort(sum(allgo(:,floor(size(allgo,2)/2):end),2) - sum(allgo(:,1:floor(size(allgo,2)/2)),2));
        [~,sorti_nogo] = sort(sum(allnogo(:,floor(size(allnogo,2)/2):end),2) - sum(allnogo(:,1:floor(size(allnogo,2)/2)),2));
        
        f = figure;
        f.Units = 'normalized';f.Position = [0.25 0.25 0.5 0.4];
        if ~indv_sorting
            f.Name = sprintf('%s sorted based on %s responses',targetcellname,sortingtype);
            s1=subplot(1,2,1);imagesc(timeoffset + xpoints,[],allgo(eval(sprintf('sorti_%s',sortingtype)),:));
        else
            f.Name = sprintf('%s sorted separately for different stimuli',targetcellname);
            s1=subplot(1,2,1);imagesc(timeoffset + xpoints,[],allgo(eval(sprintf('sorti_%s','go')),:));
        end
        greymap=colormap('gray');
        colormap(flipud(greymap));
        set(gca,'CLim',[0 clm]);
        xlim([-100 650]);
        s1.Title.String = sprintf('%s - go',targetcellname);
        s1.Title.FontWeight = 'normal';
        s1.Title.FontSize = 9;
        s1.FontSize=8;
        s1.XLabel.String = 'Time (ms)';
        s1.YLabel.String = 'Cells';
        
        
        if ~indv_sorting
            s2=subplot(1,2,2);imagesc(timeoffset + xpoints,[],allnogo(eval(sprintf('sorti_%s',sortingtype)),:));
        else
            s2=subplot(1,2,2);imagesc(timeoffset + xpoints,[],allnogo(eval(sprintf('sorti_%s','nogo')),:));
        end
        greymap=colormap('gray');
        colormap(flipud(greymap));
        set(gca,'CLim',[0 clm]);
        xlim([-100 650]);
        s2.Title.String = sprintf('%s - nogo',targetcellname);
        s2.Title.FontWeight = 'normal';
        s2.Title.FontSize = 9;
        s2.FontSize=8;
        s2.XLabel.String = 'Time (ms)';
        s2.YLabel.String = 'Cells';
    end
    %% selectivity zscore
    %% difference plots - Version2
    if plotstodo.sel       
        cellrand=nan(1000,size(targetcell{1}.nbs.go,2));% 1000* time
        %
        for rrand = 1:1000
            onerep_gonodif_percell_overtime = nan(numel(targetcell),size(targetcell{1}.nbs.go,2));
            for i= 1:numel(targetcell)
                
                target_pool = [targetcell{i}.nbs.go;targetcell{i}.nbs.nogo];
                
                selectcell=zeros(1,size(targetcell{i}.nbs.go,1)+size(targetcell{i}.nbs.nogo,1))';
                selectcell(randi(length(selectcell),[1,size(targetcell{i}.nbs.nogo,1)]))=1;
                % gonogodif percell
                if selindex
                    onerep_gonodif_percell_overtime(i,:) = ((nanmean(target_pool((find(~selectcell)),:))-nanmean(target_pool((find(selectcell)),:)))./ ...
                        (nanmean(target_pool((find(~selectcell)),:))+nanmean(target_pool((find(selectcell)),:))) )...
                        ;
                else
                    onerep_gonodif_percell_overtime(i,:) = (nanmean(target_pool((find(~selectcell)),:))-nanmean(target_pool((find(selectcell)),:)))...
                        ;
                end
            end
            cellrand(rrand,:) = nanmean(zscore(onerep_gonodif_percell_overtime')',1) ;
        end
        %%%
        if selindex
            gonogodifpercell = (allgo-allnogo)./(allgo+allnogo);
        else
            allgon=(allgo-(cellfun(@(x) x.Avfr,targetcell))')./(cellfun(@(x) x.Avfr,targetcell))';
            allnogon=(allnogo-(cellfun(@(x) x.Avfr,targetcell))')./(cellfun(@(x) x.Avfr,targetcell))';
            gonogodifpercell = (allgon-allnogon)./(allgon+allnogon);
            
        end
        figure;plot(timeoffset + xpoints,(nanmean(zscore(gonogodifpercell')',1)),'Color','b','LineWidth',1.6);
        patch_xpoints=[timeoffset + xpoints fliplr(timeoffset + xpoints)];
        patch_ypoints=[prctile(cellrand,2.5),fliplr(prctile(cellrand,97.5))];
        hold on;patch(patch_xpoints,patch_ypoints,'b','FaceAlpha',0.2,'EdgeAlpha',0)
        
        xlabel('Time (ms)');
        xlim([-100 650]);ylim([-1 1]);
        axx = gca;
        hold on; patch([0 stimlengthinms stimlengthinms 0],[axx.YLim(1) axx.YLim(1) axx.YLim(2) axx.YLim(2)],[0.7 0.7 0.7],...
            'FaceAlpha',0.3,'EdgeColor','none')
        set(gcf,'Color','w')
        axx.XTick=0:100:650;
        axx.XTickLabel={'0','','200','','400','','600',''};
    end
end
