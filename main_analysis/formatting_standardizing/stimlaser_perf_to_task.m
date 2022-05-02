% For VL animals (both FF and FB), stimlaser and perfiles are used to creat a task_.. file
% with the same format as MPV17,18,... which is saved in the same
% preprocessing folder

clear all

animallist = {'VL53','VL52','VL51','VL66',...
    'VL61','VL63','VL55','VL59','VL50'};
preprocessinglist = {'2017_12_05_18_58_17','2017_12_04_15_08_47','2017_12_03_17_30_15','2018_03_20_11_53_19',...
    '2018_03_16_16_53_13','2018_04_01_12_24_38','2018_04_01_15_03_48','2018_04_01_13_24_48','2017_12_01_16_21_49'};


for animal_i=1:length(animallist) % 1
    cd('/mnt/data/Mitra/figs');
    cd(sprintf('%s/preprocessing/%s',animallist{animal_i},preprocessinglist{animal_i}))
    
    matfilename = dir('s*.mat');
    perfmatfilename = dir('p*.mat');
    load(perfmatfilename.name)
    load(matfilename.name)
    
    %%%%%%%%%%%%%%%
    PAllOn = PAllOn{1};
    rec_length = rec_length{1};
    LAllOn = LAllOn{1};
    LStepTimeOn = LStepTimeOn{1};
    LStepTimeStampOff = LStepTimeStampOff{1};
    LStepTimeStampOn = LStepTimeStampOn{1};
    path2kwd = path2kwd{1};
    path2kwe =  path2kwe{1};
    PStepTimeStampOff = PStepTimeStampOff{1};
    PStepTimeStampOn = PStepTimeStampOn{1};
    SSmessages = SSmessages{1};
    start_time = start_time{1};
    
    nogroomingind = nogroomingind{1};
    nogotrialind = nogotrialind{1};
    gotrialind = gotrialind{1};
    correctgotrialind = correctgotrialind{1};
    correctnogotrialind = correctnogotrialind{1};
    incorrectgotrialind = incorrectgotrialind{1};
    incorrectnogotrialind = incorrectnogotrialind{1};
    %%%%%%%%%%%%%%
    
    filename = 'task_remade_from_stimlaser_perf.mat';
    copyfile(matfilename.name,filename);
    save(filename,'PAllOn','rec_length','LAllOn','LStepTimeOn',...
        'LStepTimeStampOff','LStepTimeStampOn','path2kwd','path2kwe',...
        'PStepTimeStampOff','PStepTimeStampOn','SSmessages','start_time',...
        'nogroomingind','nogotrialind','gotrialind','correctgotrialind' ,...
        'correctnogotrialind','incorrectgotrialind','incorrectnogotrialind',...
        '-v7.3','-append');
    
end