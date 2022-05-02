animallist = {'MPV17','MPV18_2',...
    'VL53','VL52','VL51','VL66'};
preprocessinglist = {'2019_04_05_21_25_20','2019_04_05_15_04_39',...
    '2017_12_05_18_58_17','2017_12_04_15_08_47','2017_12_03_17_30_15','2018_03_20_11_53_19'};
exptype = {'FB','FB',...
    'FB','FB','FB','FB'};

trialtypelist = {'gotrialind','nogotrialind'};
targettrialtypelist = {'correctgotrialind','correctnogotrialind'};
trialtype_namelist = {'go','nogo'};

allgolicks=nan(1,2*Window+1);
allnogolicks=nan(1,2*Window+1);
for animal_i=1:length(animallist)
    cd('/mnt/data/Mitra/figs');
    cd(sprintf('%s/preprocessing/%s',animallist{animal_i},preprocessinglist{animal_i}))
    if animal_i<3
        matfilename = dir('task_*.mat');
    else
        matfilename = dir('s*.mat');
        perfmatfilename = dir('p*.mat');
        load(perfmatfilename.name)
    end
    load(matfilename.name)
    if animal_i>2
        PAllOn = PAllOn{1};
        %licks = licks{1};
        nogroomingind = nogroomingind{1};
        nogotrialind = nogotrialind{1};
        gotrialind = gotrialind{1};      
        correctgotrialind = correctgotrialind{1};
        correctnogotrialind = correctnogotrialind{1};
        incorrectgotrialind = incorrectgotrialind{1};
        incorrectnogotrialind = incorrectnogotrialind{1};
    end
    size(intersect(gotrialind,nogroomingind))
    allgolicks= [allgolicks ;licks(intersect(gotrialind,nogroomingind),:)];
    allnogolicks= [allnogolicks ;licks(intersect(nogotrialind,nogroomingind),:)];
    
end
%%
figure;
beh_f2=gcf;beh_f2.Name= 'excluding grooming';
beh_f2.Units = 'normalized';
beh_f2.Position = [0.3 0.4 0.5 0.4];
subplot(1,2,1);
ax1=imagesc(((-Window+1)/samplingf:(speedWindowms/1000):Window/samplingf),[],...
    [allgolicks;nan(1,size(licks,2));allnogolicks]);
colormap(ax1.Parent,flipud(gray(10)))
ax1.Parent.CLim = [0 1];
ax1.Parent.XLabel.String='time(s)';
ax1.Parent.YLabel.String='trials';
ax1.Parent.Title.String='licks aligned to stimulus onset';
