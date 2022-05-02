function [ALaserTimeStamp,LaserDelay,centers,LaserDelayBinned] = ExtractlaserdelayBeforeStim(LStepTimeOn,PDAllOn,Window,samplingf,laserbinms,laserhistbinnum)

pdtriged_Laser = LStepTimeOn(PDAllOn(:,:));
figure;plot((Window(1):Window(2))/samplingf,pdtriged_Laser');
%%% absolute laser times
middleTimeStamp=nan(1,size(PDAllOn,1));
for i=1:length(middleTimeStamp)
    if numel(find(pdtriged_Laser(i,:)))>0
    middleTimeStamp(i)=(PDAllOn(i,find(pdtriged_Laser(i,:))));
    end
end
ALaserTimeStamp= repmat(middleTimeStamp',[1,(Window(2)-Window(1))+1])+repmat(double(Window(1):Window(2)),[length(middleTimeStamp'),1]); 

%%%
LaserDelay = nan(1,size(pdtriged_Laser,1));
for i= 1: size(pdtriged_Laser,1)
    if numel(find(pdtriged_Laser(i,:)))
        LaserDelay(i) = (find(pdtriged_Laser(i,:),1) + Window(1) ) / (samplingf/1000) ; %in ms
    end
end

figure;hist(LaserDelay,laserhistbinnum)

%%%
[b,h] = histcounts(LaserDelay,laserhistbinnum); %200
[a,i]=findpeaks([ 0 0 b 0 0]);
centers=h(i - 2);


LaserDelayBinned = zeros(size(LaserDelay));
for i = 1:length(centers)    
     LaserDelayBinned(LaserDelay >= centers(i)-laserbinms & LaserDelay < centers(i)+laserbinms) = i; % 1.5 instead of 1
end
figure;hist(LaserDelayBinned,laserhistbinnum)