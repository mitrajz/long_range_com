
min_n_trials_per_delay = 20;
excludegrooming=0; 

animallist = {'VL53','VL52','VL51','VL66',...
    'VL61','VL63','VL55','VL59','VL50'};
preprocessinglist = {'2017_12_05_18_58_17','2017_12_04_15_08_47','2017_12_03_17_30_15','2018_03_20_11_53_19',...
    '2018_03_16_16_53_13','2018_04_01_12_24_38','2018_04_01_15_03_48','2018_04_01_13_24_48','2017_12_01_16_21_49'};
exptype = {'FB','FB','FB','FB',...
    'FF','FF','FF','FF','FF'};

respwindow = [100 550];
FA=nan(1,length(animallist));
Miss=nan(1,length(animallist));

for animal_i=1:length(animallist)

    cd('/mnt/data/Mitra/figs');
    cd(sprintf('%s/preprocessing/%s',animallist{animal_i},preprocessinglist{animal_i}))
    matfilename = dir('task_*.mat');
    load(matfilename.name)
    
    correctnogotrialind = [];
    correctgotrialind = [];
    incorrectnogotrialind = [];
    incorrectgotrialind = [];
    
    FA(animal_i) = 0;
    Miss(animal_i) = 0;
    if excludegrooming
        indss_nogo = intersect(nogotrialind,nogroomingind);
        indss_go = intersect(gotrialind,nogroomingind);
    else
        indss_nogo = nogotrialind;
        indss_go = gotrialind;
    end
    for i = 1:length(indss_nogo)
        firstlick = find(licks(indss_nogo(i),floor(size(licks,2)/2):end),1)/30;
        if firstlick>respwindow(1) & firstlick<respwindow(2)
            FA(animal_i)=FA(animal_i)+1;
            incorrectnogotrialind(end+1) = indss_nogo(i);
        elseif firstlick<respwindow(1)
        else
            correctnogotrialind(end+1) = indss_nogo(i);
        end
    end
    
    for i = 1:length(indss_go)
        firstlick = find(licks(indss_go(i),floor(size(licks,2)/2):end),1)/30;
        if firstlick>respwindow(1) & firstlick<respwindow(2)
            correctgotrialind(end+1) = indss_go(i);
        elseif firstlick<respwindow(1)
        else
            Miss(animal_i)=Miss(animal_i)+1;
            incorrectgotrialind(end+1) = indss_go(i);
        end
    end
    
    FA(animal_i)= 100*FA(animal_i)/length(indss_nogo)
    Miss(animal_i)=100*Miss(animal_i)/length(indss_go)
    
    save('rev_perf.mat','correctnogotrialind','correctgotrialind','incorrectnogotrialind','incorrectgotrialind');
end

