
x=[]
load('cells_FF_pl20_an150_lw20_exG1_onlyC1_onlyS0_plstyle1.mat')
for i=unique(cellfun(@(x) x.simulcode,V1cells))
    
    x=[x; V1cells{find(cellfun(@(x) x.simulcode,V1cells) ==i,1)}.smb_centers.go];
end

load('cells_FB_pl20_an150_lw20_exG1_onlyC1_onlyS0_plstyle1.mat')
for i=unique(cellfun(@(x) x.simulcode,V1cells))
    
    x=[x; V1cells{find(cellfun(@(x) x.simulcode,V1cells) ==i,1)}.smb_centers.go];
end

round(mean(x,1)) % mean across animals
(round(std(x,1))) % std across animals

%% jitter within animal -- jitter

animallist ={'VL61','VL63','VL55','VL59',...
    'MPV33','MPV31','MPV34_2',...
    'MPV17','MPV18_2',...
    'VL53','VL52','VL51','VL66','MPV35_2'};
preprocessinglist = {'2018_03_16_16_53_13','2018_04_01_12_24_38','2018_04_01_15_03_48','2018_04_01_13_24_48',...
    '2020_03_02_19_35_12','2020_03_02_19_58_50','2020_03_02_20_32_35',...
    '2020_03_02_13_10_37','2020_03_02_13_27_45',...
     '2017_12_05_18_58_17','2017_12_04_15_08_47','2017_12_03_17_30_15','2018_03_20_11_53_19','2020_03_02_14_31_53',...
     };

exptype = {'FF','FF','FF','FF',...
    'FF','FF','FF',...
    'FB','FB',...
    'FB','FB','FB','FB','FB'};

range_ms = nan(length(animallist),8);
std_ms = nan(length(animallist),8);
n_eachlag_trials = nan(length(animallist),8);
n_nolaser_trials = nan(length(animallist),1);

for animal_i=1:length(animallist) 
    cd('/mnt/data/Mitra/figs');
    cd(sprintf('%s/preprocessing/%s',animallist{animal_i},preprocessinglist{animal_i}))
    matfilename = dir('task_*.mat');
    load(matfilename.name)
    
    lagcount = 0;
    for i = unique(LaserDelayBinned(find(~isnan(LaserDelayBinned))))        
        if numel(find(LaserDelayBinned == i))>40 % this is full trials, 
                                                 % so, a full session should have around 80
            lagcount = lagcount+1;
            std_ms(animal_i,lagcount) = std(LaserDelay(find(LaserDelayBinned == i)));
            range_ms(animal_i,lagcount) = range(LaserDelay(find(LaserDelayBinned == i)));
            n_eachlag_trials(animal_i,lagcount) = numel(LaserDelay(find(LaserDelayBinned == i)));
            
        end
    end
    n_nolaser_trials(animal_i) = numel(LaserDelay(find(isnan(LaserDelayBinned))));
end

% max std
max(max(std_ms))

% average number of trials per animal, before any exclusion
mean(reshape(n_eachlag_trials,1,[]))
std(reshape(n_eachlag_trials,1,[]))
% average numer of no laser trials
mean(n_nolaser_trials)
std(n_nolaser_trials)