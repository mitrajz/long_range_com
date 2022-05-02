function [fouthandle,LaserDelay,centers,LaserDelayBinned] = Extractlaserdelaynan(LStepTimeOn,PDAllOn,Window,samplingf,laserbinms,laserhistbinnum)
% LaserDelay is the true laser delay of each trial. For any trial with a
% laser it is a number in ms and for non-laser trials it is nan. The only
% correct way of picking no laser trials is find(isnan(LaserDelay))
%
% centers is the center in ms of the found peaks. There is no guarantee on
% the number of laser trials in each center group. 
%
% LaserDelayBinned
%
% when using method2 (default), laserbinms could be nan, which means there
% is no limit on center assignment, otherwise, it can be set to a number in
% ms which determines the maximum jitter in each center and discards the
% trials whose delay is not within within this range from any center
%
%
% This version is improved a lot compared to Extractlaserdelay
%
pdtriged_Laser = LStepTimeOn(PDAllOn(:,:));
fouthandle = figure('Units','normalized','Position',[0.2 0.1 0.7 0.7]);
s1 = subplot(3,2,1); plot((-Window:Window)/samplingf,pdtriged_Laser');
s1.Title.String = 'LaserDelay(aligned to pd) for all trials';

LaserDelay = nan(1,size(pdtriged_Laser,1));
for i= 1: size(pdtriged_Laser,1)
    temp = find(pdtriged_Laser(i,:));
   % temp = (temp(temp>=(size(PDAllOn,2)-1)/2)) ;
     temp = (temp(temp>=(size(PDAllOn,2)-1)/2 - 20*30)) ; % to get up to 20
    % ms before 
    if numel(temp)
        
        LaserDelay(i) =(temp - (size(PDAllOn,2)-1)/2) / (samplingf/1000) ; %in ms
    end
end

s2 = subplot(3,2,2); hist(LaserDelay,laserhistbinnum)
s2.Title.String = 'histogram of laserDelays with the given bin';

%%%
[b,h] = histcounts(LaserDelay,laserhistbinnum); %200
[a,i]=findpeaks([ 0 0 b 0 0]);
centers=h(i - 2);

LaserDelayBinned = nan(size(LaserDelay));
method = 2;

if method == 1
  
    for i = 1:length(centers)
        LaserDelayBinned(LaserDelay >= centers(i)-laserbinms & LaserDelay < centers(i)+laserbinms) = i;
    end
elseif method == 2
    % This overcomes the problem of the first method by assigning each
    % laser delay to the center closest to it. 
    for i = 1:length(LaserDelay)
        if ~isnan(LaserDelay(i))
            [min_ms,closest_center_ind] = min(abs(LaserDelay(i) - centers));
            if ~isnan (laserbinms)
                if min_ms < laserbinms
                    LaserDelayBinned(i)  = closest_center_ind;
                end
                % now, centers(closest_center_ind) will give the binned delay
                % in ms of the trial
            else
                LaserDelayBinned(i)  = closest_center_ind;                
            end
        end
    end
    
end
s3 = subplot(3,2,3); hist(LaserDelayBinned,length(centers));
s3.Title.String = 'LaserDelayBinned';

%calculate jitter within each center
jitter = nan(length(centers));
for i=1:length(centers)
    if numel(find(LaserDelayBinned == i)) > 0
        jitter(i) = max(LaserDelay(find(LaserDelayBinned == i))) - min(LaserDelay(find(LaserDelayBinned == i)));
    end
end
s4 = subplot(3,2,4); plot(centers , jitter);
s4.Title.String = 'jitter for each center in ms';
%% checks:
s5=subplot(3,2,[5,6]);
text(0.1,0.5,sprintf('\nnumber of nolaser trials = %d\nand number of laser trials = %d\n%dcenetrs found, of which %d are non empty\n%d laser trials were not binned\nexpected to find around %d centers\n',...
    numel(find(isnan(LaserDelay))),numel(find(~isnan(LaserDelay))),...
    length(centers),length(unique(LaserDelayBinned(find(~isnan(LaserDelayBinned))))),...
    numel(find(isnan(LaserDelayBinned))) - numel(find(isnan(LaserDelay))),...
    round((centers(end) - centers(1)) / abs(mode(diff(LaserDelay)))) ));
set( s5, 'visible', 'off')

fprintf('\nnumber of nolaser trials = %d\nand number of laser trials = %d\n',numel(find(isnan(LaserDelay))),numel(find(~isnan(LaserDelay))));
fprintf('%d cenetrs found, of which %d are non empty\n',length(centers),length(unique(LaserDelayBinned(find(~isnan(LaserDelayBinned))))));
fprintf('%d laser trials were not binned\n',numel(find(isnan(LaserDelayBinned))) - numel(find(isnan(LaserDelay))));
fprintf('expected to find around %d centers\n', round((centers(end) - centers(1)) / abs(mode(diff(LaserDelay)))));

