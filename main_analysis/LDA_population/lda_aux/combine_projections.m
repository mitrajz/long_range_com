function combine_projections(proj_bs_go,proj_ls_go,proj_dif_go,proj_animalcount_bs_go,proj_animalcount_ls_go,...
    proj_bs_nogo,proj_ls_nogo,proj_dif_nogo,proj_animalcount_bs_nogo,proj_animalcount_ls_nogo)
%%
coef = 1;

param = 2;
method = 2; % 
subtractprestim = 1;
% subtract bias
% this is individual
for i=1:8
    if method == 1
        proj_bs_go{i} = (proj_bs_go{i}-nanmean(proj_bs_go{i}(:,1:size(proj_bs_go{i},2)/param),2));
        proj_ls_go{i} = (proj_ls_go{i}-nanmean(proj_ls_go{i}(:,1:size(proj_ls_go{i},2)/param),2));
        proj_bs_nogo{i} = (proj_bs_nogo{i}-nanmean(proj_bs_nogo{i}(:,1:size(proj_bs_nogo{i},2)/param),2));
        proj_ls_nogo{i} = (proj_ls_nogo{i}-nanmean(proj_ls_nogo{i}(:,1:size(proj_ls_nogo{i},2)/param),2));
    elseif method == 2
        if subtractprestim % test, for each trial
            proj_bs_go{i} = (proj_bs_go{i}-nanmean(proj_bs_go{i}(:,1:size(proj_bs_go{i},2)/param),2));
            proj_ls_go{i} = (proj_ls_go{i}-nanmean(proj_ls_go{i}(:,1:size(proj_ls_go{i},2)/param),2));
            proj_bs_nogo{i} = (proj_bs_nogo{i}-nanmean(proj_bs_nogo{i}(:,1:size(proj_bs_nogo{i},2)/param),2));
            proj_ls_nogo{i} = (proj_ls_nogo{i}-nanmean(proj_ls_nogo{i}(:,1:size(proj_ls_nogo{i},2)/param),2));
        end
        for animali = unique(proj_animalcount_bs_go{i})'
            ind_bs = find(proj_animalcount_bs_go{i} == animali);
            ind_ls = find(proj_animalcount_ls_go{i} == animali);
            
  
            bias = nanmean(nanmean([proj_bs_go{i}(ind_bs,:);proj_ls_go{i}(ind_ls,:)]));
            proj_bs_go{i}(ind_bs,:) = (proj_bs_go{i}(ind_bs,:) -bias);
            proj_ls_go{i}(ind_ls,:) = (proj_ls_go{i}(ind_ls,:) -bias);
            
        end
        for animali = unique(proj_animalcount_bs_nogo{i})'
            ind_bs = find(proj_animalcount_bs_nogo{i} == animali);
            ind_ls = find(proj_animalcount_ls_nogo{i} == animali);
            
            
   
            bias = nanmean(nanmean([proj_bs_nogo{i}(ind_bs,:);proj_ls_nogo{i}(ind_ls,:)]));
            proj_bs_nogo{i}(ind_bs,:) = (proj_bs_nogo{i}(ind_bs,:) -bias);
            proj_ls_nogo{i}(ind_ls,:) = (proj_ls_nogo{i}(ind_ls,:) -bias);
        end
        
    end
    
    
end


smoothspan = 5;

figure;
for i=1:8
    subplot(8,1,i);
    shadedErrorBar([],smooth(nanmean(proj_bs_go{i},1),smoothspan),...
        smooth(coef*nanstd(proj_bs_go{i})/sqrt(size(proj_bs_go{i},1)),smoothspan),'k',1);
    hold on;
    shadedErrorBar([],smooth(nanmean(proj_ls_go{i},1),smoothspan),...
        smooth(coef*nanstd(proj_ls_go{i})/sqrt(size(proj_ls_go{i},1)),smoothspan),'b',1)
    
end


figure;
for i=1:8
    subplot(8,1,i);
    shadedErrorBar([],smooth(nanmean(proj_bs_nogo{i},1),smoothspan),...
        smooth(coef*nanstd(proj_bs_nogo{i})/sqrt(size(proj_bs_nogo{i},1)),smoothspan),'k',1);
    hold on;
    shadedErrorBar([],smooth(nanmean(proj_ls_nogo{i},1),smoothspan),...
        smooth(coef*nanstd(proj_ls_nogo{i})/sqrt(size(proj_ls_nogo{i},1)),smoothspan),'b',1)
    
end


