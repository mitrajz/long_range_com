% This script is to be run after csd and depth assignment for individual animals. 
% For each shank, csd and lfp of all animals are plotted together to check
% alignment. If for an animal, depth assignment wa not consistent with the
% others, run the 2nd section to correct and overwite it.

animallist = {'MPV17','MPV18_2',...
    'VL53','VL52','VL51','VL66'};%,...
%'VL61','VL63','VL55','VL59','VL50'};
preprocessinglist = {'2019_04_05_21_25_20','2019_04_05_15_04_39',...
    '2017_12_05_18_58_17','2017_12_04_15_08_47','2017_12_03_17_30_15','2018_03_20_11_53_19',...
    '2018_03_16_16_53_13','2018_04_01_12_24_38','2018_04_01_15_03_48','2018_04_01_13_24_48','2017_12_01_16_21_49'};

% animallist = {'UMPV2','UMPV4','UMPV7','UMPV8'};
% preprocessinglist = {'2019_09_13_17_56_13','2019_09_13_18_12_17',...
%     '2019_09_13_18_35_50','2019_09_13_19_01_09'};

% checking and plotting 1 shank at a time for all animals;
shankname = 'LM_1';
allf = figure;
allf.Name = shankname;
for animal_i=1:length(animallist) % 1
    allsub_csd=subplot(length(animallist),2,2*animal_i-1,'Parent',allf);
    allsub_lfp=subplot(length(animallist),2,2*animal_i,'Parent',allf);
    cd('/mnt/data/Mitra/figs');
    cd(sprintf('%s/preprocessing/%s',animallist{animal_i},preprocessinglist{animal_i}))
    load('shankinfo.mat')
    uiopen(['lfp_csd_depth_',shankname,'.fig'],1);
    shanknum=find(strcmp(order,shankname));
    savedf = gcf;
    smoothlfp = savedf.Children(3);
    smoothcsd = savedf.Children(2);
    %%%% fill in csd
    copyobj(smoothcsd.Children,allsub_csd)
    allsub_csd.YDir='reverse';
    allsub_csd.YTick=refdepth{shanknum};%1:4:numel(find(strcmp(animalshankorder,order{shanknum})))-2;
    allsub_csd.YTickLabel= '0';%shankdepthmap{shanknum}(1:4:end-2);
    hold(allsub_csd,'on');line(allsub_csd,[-50 200] , [refdepth{shanknum}, refdepth{shanknum}],'Color','k')
    hold(allsub_csd,'on');line(allsub_csd,[-50 200] , [refdepth{shanknum}-4, refdepth{shanknum}-4],'Color','k') % 100 microns above: roughly top of L4
    allsub_csd.XLim = [-50 200];
    allsub_csd.YLim = [0 32];
    %allsub_csd.CLim = [-40 40];
    %%%% fill in lfp
    copyobj(smoothlfp.Children,allsub_lfp)
    allsub_lfp.YDir='reverse';
    allsub_lfp.YTick=refdepth{shanknum};%1:4:numel(find(strcmp(animalshankorder,order{shanknum})))-2;
    allsub_lfp.YTickLabel= '0';%shankdepthmap{shanknum}(1:4:end-2);
    hold(allsub_lfp,'on');line(allsub_lfp,[-50 200] , [refdepth{shanknum}, refdepth{shanknum}],'Color','k')
    hold(allsub_lfp,'on');line(allsub_lfp,[-50 200] , [refdepth{shanknum}-4, refdepth{shanknum}-4],'Color','k') % 100 microns above: roughly top of L4
    allsub_lfp.XLim = [-50 200];
    allsub_lfp.YLim = [0 32];
    % allsub_lfp.CLim = [-600 600];
end
%% redo if needed
% if for 1 animal, a certain shank had incorrect depth aignment, run this
% part. set animal_i (according to list above), shankname to redo, and set
% savevars to 1 if you want to overwrite depth info

if false
    animal_i = 2;
    shank_to_redo_name = 'LM_1';
    smoothinglength = 5;
    savevars = 0;
    nchannels = 128;
    
    cd('/mnt/data/Mitra/figs');
    cd(sprintf('%s/preprocessing/%s',animallist{animal_i},preprocessinglist{animal_i}))
    load('shankinfo.mat')
    load('stimalignedlfp.mat', 'data_lfp_bs')
    load('stimalignedlfp.mat', 'lfptimewindow')
    animalshankorder = cell(1,nchannels);
    for i=1:nchannels
        if i<=nchannels/4
            animalshankorder{i} = order{1};
        elseif i<=nchannels/2
            animalshankorder{i} = order{2};
        elseif i<=3*nchannels/4
            animalshankorder{i} = order{3};
        else
            animalshankorder{i} = order{4};
        end
    end
    % only redoes for 1 shank. Not thouroughly tested, but should work.
    shanknum=find(strcmp(order,shank_to_redo_name));
    [shankdepthmap(shanknum),refdepth(shanknum),lfpdepth_f(shanknum)] = alignlfpandgetdepth(order(shanknum),data_lfp_bs,animalshankorder,lfptimewindow,smoothinglength);
   
    % after you open all shanks
    if savevars
        cd('/mnt/data/Mitra/figs');
        cd(sprintf('%s/preprocessing/%s',animallist{animal_i},preprocessinglist{animal_i}))
        % save depth to shankinfo.mat
        save('shankinfo.mat','refdepth','shankdepthmap','-append');
        % save lfpdepth figures
        %for shanknum = 1:numel(order)
        shanknum = find(strcmp(order,shank_to_redo_name));
        saveas(lfpdepth_f{shanknum},sprintf('lfp_csd_depth_%s.fig',order{shanknum}))
       % end
        close all
    end
end
