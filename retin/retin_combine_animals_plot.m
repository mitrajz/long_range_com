% make sure to run retin_pipeline.m for all animals in the list before running this script. 
%%
clear all
%%
addpath(genpath('/mnt/data/Mitra/cache/codes/matlabkwikreader/'))
addpath(genpath('/mnt/data/Mitra/cache/codes/'))
%%
% animallist = {'MPV17','MPV18_2',...
%     'VL53','VL52','VL51','VL66',...
%     'VL61','VL63','VL55','VL59','VL50',...
%     'MPV32_2','MPV33','MPV31','MPV34_2',...
%     'MPV36_2','MPV30','MPV35_2'};
% preprocessinglist = {'2019_04_05_21_25_20','2019_04_05_15_04_39',...
%     '2017_12_05_18_58_17','2017_12_04_15_08_47','2017_12_03_17_30_15','2018_03_20_11_53_19',...
%     '2018_03_16_16_53_13','2018_04_01_12_24_38','2018_04_01_15_03_48','2018_04_01_13_24_48','2017_12_01_16_21_49',...
%     '2020_01_21_18_13_31','2020_01_21_17_12_32','2020_01_22_13_47_56','2020_01_22_14_28_52',...
%     '2020_01_22_11_38_57','2020_01_22_14_51_15','2020_01_28_20_29_59'};

animallist = {'MPV17','MPV18_2',...
    'VL53','VL52','VL51','VL66',...
    'VL61','VL63','VL55','VL59',...
    'MPV33','MPV31','MPV34_2',...
    'MPV36_2','MPV30','MPV35_2'};
preprocessinglist = {'2019_04_05_21_25_20','2019_04_05_15_04_39',...
    '2017_12_05_18_58_17','2017_12_04_15_08_47','2017_12_03_17_30_15','2018_03_20_11_53_19',...
    '2018_03_16_16_53_13','2018_04_01_12_24_38','2018_04_01_15_03_48','2018_04_01_13_24_48',...
    '2020_01_21_17_12_32','2020_01_22_13_47_56','2020_01_22_14_28_52',...
    '2020_01_22_11_38_57','2020_01_22_14_51_15','2020_01_28_20_29_59'};

Normalize = 0; %normalize each units map to its max or not, just for plotting

all_V1s = cell(0);
all_LMs = cell(0);
for animal_i=1:length(animallist)
    
    animallist{animal_i}
    
    cd('/mnt/data/Mitra/figs');
    cd(sprintf('%s/preprocessing/%s',animallist{animal_i},preprocessinglist{animal_i}))
    
    if ~exist('retinmaps.mat','file')
        disp('no retin map');
    else        
        load('retinmaps.mat');
        all_V1s{end+1} = retin_V1{1}.retinmap;
        all_V1s{end+1} = retin_V1{2}.retinmap;
        all_LMs{end+1} = retin_LM{1}.retinmap;
        all_LMs{end+1} = retin_LM{2}.retinmap;        
    end
end

f_m=figure;
mergeshanksandplot(all_V1s,'taskmatched',subplot(2,2,1),subplot(2,2,3),'V1',Normalize)
mergeshanksandplot(all_LMs,'taskmatched',subplot(2,2,2),subplot(2,2,4),'LM',Normalize)
f_m.Name = [animallist{animal_i},'-','taskmatched'];

f_g=figure;
mergeshanksandplot(all_V1s,'KSgood',subplot(2,2,1),subplot(2,2,3),'V1',Normalize)
mergeshanksandplot(all_LMs,'KSgood',subplot(2,2,2),subplot(2,2,4),'LM',Normalize)
f_g.Name = [animallist{animal_i},'-','KSgood'];

f_all=figure;
mergeshanksandplot(all_V1s,'all',subplot(2,2,1),subplot(2,2,3),'V1',Normalize)
mergeshanksandplot(all_LMs,'all',subplot(2,2,2),subplot(2,2,4),'LM',Normalize)
f_all.Name = [animallist{animal_i},'-','all'];
