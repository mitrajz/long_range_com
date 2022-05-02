x=1;
min_n_trials_per_delay = 20;

animallist = {'VL53','VL52','VL51','VL66','VL58',...
    'VL61','VL63','VL55','VL59','VL50'};
preprocessinglist = {'2017_12_05_18_58_17','2017_12_04_15_08_47','2017_12_03_17_30_15','2018_03_20_11_53_19','2018_03_26_15_29_09',...
    '2018_03_16_16_53_13','2018_04_01_12_24_38','2018_04_01_15_03_48','2018_04_01_13_24_48','2017_12_01_16_21_49'};
exptype = {'FB','FB','FB','FB','FB',...
    'FF','FF','FF','FF','FF'};
indss_go = cell(1,length(animallist));
indss_nogo = cell(1,length(animallist));
for animal_i=1:length(animallist)
    %cd('/mnt/data/Mitra/figs/P2_L/LED/FB');
    cd('/mnt/data/Mitra/figs');
    cd(sprintf('%s/preprocessing/%s',animallist{animal_i},preprocessinglist{animal_i}))
    matfilename = dir('*.mat');
    load(matfilename.name)
    
    
    %%%%%%%%%
    % indices of baseline trials
    
    %indss0=intersect(intersect(find(isnan(LaserDelay)),eval(trialtype)),nogroomingind{x});
    indss_go{animal_i}=nan(1,length(centers));
    indss_nogo{animal_i}=nan(1,length(centers));
    for ldelay =1:length(centers)
        indss_g= intersect(find(LaserDelayBinned == (ldelay)),gotrialind{x});
        indss_ng= intersect(find(LaserDelayBinned == (ldelay)),nogotrialind{x});
        if length(indss_g) > min_n_trials_per_delay
            indss_go{animal_i}(1,ldelay)=length(intersect(indss_g,nogroomingind{x}));
        end
        if length(indss_ng) > min_n_trials_per_delay
            indss_nogo{animal_i}(1,ldelay)=length(intersect(indss_ng,nogroomingind{x}));
        end
    end
    
    indss_go{animal_i} = indss_go{animal_i}(find(~isnan(indss_go{animal_i})));
    indss_nogo{animal_i} = indss_nogo{animal_i}(find(~isnan(indss_nogo{animal_i})));
    
end

indss_nogo_mat = cell2mat(indss_nogo');
indss_go_mat = cell2mat(indss_go');