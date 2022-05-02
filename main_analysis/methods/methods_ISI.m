
% only trained animals here 

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

ISI_sec = cell(1,length(animallist));
ISI_cap = nan(1,length(animallist));


for animal_i=1:length(animallist) 
    cd('/mnt/data/Mitra/figs');
    cd(sprintf('%s/preprocessing/%s',animallist{animal_i},preprocessinglist{animal_i}))
    matfilename = dir('task_*.mat');
    load(matfilename.name)
    
    
    ISI_sec{animal_i} = diff(PStepTimeStampOn)/30000;
    ISI_cap(animal_i) = max(diff(PStepTimeStampOn)/30000);
    
end

ISI_cap
cellfun(@(x) nanmean(x), ISI_sec) % the mean is post capping

%figure; hist(ISI_sec{12},50)
