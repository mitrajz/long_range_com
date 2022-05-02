rootdir = '/mnt/data/Mitra/figs';

a = dir(fullfile(rootdir,'/**/*.fig'));
for i=1:length(a)
    i
    cd(a(i).folder)
    
    uiopen(fullfile(a(i).folder,a(i).name),1)
    f=gcf;
    set(f,'Color','w')
  
    saveas(f,[a(i).name(1:end-4),'.svg'])    
    close(f);
    pause(1)
end




