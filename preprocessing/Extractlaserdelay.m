function [LaserDelay,centers,LaserDelayBinned] = Extractlaserdelay(LStepTimeOn,PDAllOn,Window,samplingf,laserbinms,laserhistbinnum)
if length(Window) == 1
    Window=repmat(Window,[1,2]).*[-1, 1];
end

pdtriged_Laser = LStepTimeOn(PDAllOn(:,:));
figure;plot((Window(1):Window(2))/samplingf,pdtriged_Laser');

LaserDelay = nan(1,size(pdtriged_Laser,1));
for i= 1: size(pdtriged_Laser,1)
    if numel(find(pdtriged_Laser(i,:)))
        LaserDelay(i) = (find(pdtriged_Laser(i,:),1) - (size(PDAllOn,2)-1)/2) / (samplingf/1000) ; %in ms
    end
end

figure;hist(LaserDelay,laserhistbinnum)

%%%
[b,h] = histcounts(LaserDelay,laserhistbinnum);
[a,i]=findpeaks([ 0 0 b 0 0]);
centers=h(i - 2);


LaserDelayBinned = zeros(size(LaserDelay));
for i = 1:length(centers)    
     LaserDelayBinned(LaserDelay >= centers(i)-laserbinms & LaserDelay < centers(i)+laserbinms) = i; % 1.5 instead of 1
end
figure;hist(LaserDelayBinned,100)