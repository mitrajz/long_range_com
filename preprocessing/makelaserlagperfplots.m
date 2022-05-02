


gotrialind_w = gotrialind;
nogotrialind_w = nogotrialind;
correctgotrialind_w = correctgotrialind;
incorrectgotrialind_w = incorrectgotrialind;
correctnogotrialind_w = correctnogotrialind;
incorrectnogotrialind_w = incorrectnogotrialind;


if onlystableperiod
    gotrialind_w = intersect(gotrialind_w,stabletrialind);
    nogotrialind_w = intersect(nogotrialind_w,stabletrialind);
    correctgotrialind_w = intersect(correctgotrialind_w,stabletrialind);
    incorrectgotrialind_w = intersect(incorrectgotrialind_w,stabletrialind);
    correctnogotrialind_w = intersect(correctnogotrialind_w,stabletrialind);
    incorrectnogotrialind_w = intersect(incorrectnogotrialind_w,stabletrialind);
end
if onlynogrooming
    gotrialind_w = intersect(gotrialind_w,nogroomingind);
    nogotrialind_w = intersect(nogotrialind_w,nogroomingind);
    correctgotrialind_w = intersect(correctgotrialind_w,nogroomingind);
    incorrectgotrialind_w = intersect(incorrectgotrialind_w,nogroomingind);
    correctnogotrialind_w = intersect(correctnogotrialind_w,nogroomingind);
    incorrectnogotrialind_w = intersect(incorrectnogotrialind_w,nogroomingind);
end

%%%%%%%%%%% go 
s1=subplot(6,2,[1 3]);
s1.YAxis.Label.FontWeight = 'normal';
s1.YAxis.Label.FontSize = 9;
s1.YAxis.Label.String = 'Lick latency correct go trials'; 
s1.XAxis.Label.String = 'time (ms)';
s1.XAxis.Label.FontSize = 5;

hitrate = nan;
missrate = nan;
xpoints = nan ;
for ldelay =1:length(centers)
    indss= intersect(find(LaserDelayBinned == (ldelay)),gotrialind_w);  % only in go trials
     %indss= find(LaserDelayBinned == (ldelay));
    if length(indss) > min_n_trials_per_delay
        
        hitrate(end+1) = length(intersect(correctgotrialind_w,indss))/length(indss);
        missrate(end+1) = length(intersect(incorrectgotrialind_w,indss))/length(indss);
        xpoints(end+1) = centers(ldelay);
        %  length(indss)
        indssfirstlicks=firstlicksample(indss)/30;
        
        indssfirstlicks(find(indssfirstlicks<respwindowms(1))) =nan;
        indssfirstlicks(find(indssfirstlicks>=respwindowms(2))) =nan;
        if normalizeforcorrect
            indssfirstlicks(find(isnan(indssfirstlicks))) = [];
        end
        %[centers(ldelay) numel(indssfirstlicks)]
        hold on; plot(respwindowms(1):cumsumstep:respwindowms(2)-cumsumstep,cumsum(histcounts(indssfirstlicks,respwindowms(1):cumsumstep:respwindowms(2),'Normalization','count'))/numel(indssfirstlicks),'b')
        xlim([0 1000])
        % if ~isnan(LaserDelayBinned(ldelay))
        if centers(ldelay)> 50 & centers(ldelay)< 80
            %centers(ldelay)
            hold on;plot(respwindowms(1):cumsumstep:respwindowms(2)-cumsumstep ,cumsum(histcounts(indssfirstlicks,respwindowms(1):cumsumstep:respwindowms(2),'Normalization','count'))/numel(indssfirstlicks),'k','LineWidth',2)
            xlim([0 1000])
        end
        if centers(ldelay)<0 & centers(ldelay)>-20
            %centers(ldelay)
            hold on;plot(respwindowms(1):cumsumstep:respwindowms(2)-cumsumstep ,cumsum(histcounts(indssfirstlicks,respwindowms(1):cumsumstep:respwindowms(2),'Normalization','count'))/numel(indssfirstlicks),'r','LineWidth',1)
            xlim([0 1000])
        end
        % end
    end
end
indss0= intersect(find(isnan(LaserDelayBinned)),gotrialind_w);
length(indss0);
hitrate(1) = length(intersect(correctgotrialind_w,indss0))/length(indss0);
missrate(1) = length(intersect(incorrectgotrialind_w,indss0))/length(indss0);
indssfirstlicks=firstlicksample(indss0)/30;

indssfirstlicks(find(indssfirstlicks<respwindowms(1))) =nan;
indssfirstlicks(find(indssfirstlicks>=respwindowms(2))) =nan;
if normalizeforcorrect
    indssfirstlicks(find(isnan(indssfirstlicks))) = [];
end
hold on;plot(respwindowms(1):cumsumstep:respwindowms(2)-cumsumstep ,cumsum(histcounts(indssfirstlicks,respwindowms(1):cumsumstep:respwindowms(2),'Normalization','count'))/numel(indssfirstlicks),'g','LineWidth',1)
xlim([0 1000])

subplot(6,2,2);p1=plot(xpoints,hitrate,'+-','Color',[0 1 0]);
p1.Parent.YAxis.Label.String = 'hit rate';
p1.Parent.YAxis.Label.FontWeight='normal';
p1.Parent.YAxis.Label.FontSize=7;
p1.Parent.XLim = [-50 500];
line(p1.Parent.XLim,[hitrate(1) hitrate(1)],'Color',[0.8 0.8 0.8]);

subplot(6,2,4);p2=plot(xpoints,missrate,'+-','Color',[0.2 0.5 0.2]);
p2.Parent.YAxis.Label.String = 'missrate';
p2.Parent.YAxis.Label.FontWeight='normal';
p2.Parent.YAxis.Label.FontSize=7;
p2.Parent.XLim = [-50 500];
line(p2.Parent.XLim,[missrate(1) missrate(1)],'Color',[0.8 0.8 0.8]);


%%%%%%%%%%%%%%%% nogo
crrate = nan;
farate = nan;
s2=subplot(6,2,[5 7]);
s2.YAxis.Label.FontWeight = 'normal';
s2.YAxis.Label.FontSize = 9;
s2.YAxis.Label.String = 'Lick latency incorrect nogo trials'; 
s2.XAxis.Label.String = 'time (ms)';
s2.XAxis.Label.FontSize = 5;

for ldelay =1:length(centers)
    indss= intersect(find(LaserDelayBinned == (ldelay)),nogotrialind_w);  % only in go trials
    % indss= find(LaserDelayBinned == (ldelay));
    if length(indss) > min_n_trials_per_delay
        
        crrate(end+1) = length(intersect(correctnogotrialind_w,indss))/length(indss);
        farate(end+1) = length(intersect(incorrectnogotrialind_w,indss))/length(indss);
        %  length(indss)
        indssfirstlicks=firstlicksample(indss)/30;
        
        indssfirstlicks(find(indssfirstlicks<respwindowms(1))) =nan;
        indssfirstlicks(find(indssfirstlicks>=respwindowms(2))) =nan;
        if normalizeforcorrect
            indssfirstlicks(find(isnan(indssfirstlicks))) = [];
        end
        hold on; plot(respwindowms(1):cumsumstep:respwindowms(2)-cumsumstep,cumsum(histcounts(indssfirstlicks,respwindowms(1):cumsumstep:respwindowms(2),'Normalization','count'))/numel(indssfirstlicks),'b')
        xlim([0 1000])
        % if ~isnan(LaserDelayBinned(ldelay))
        if centers(ldelay)> 50 & centers(ldelay)< 60
            %centers(ldelay)
            hold on;plot(respwindowms(1):cumsumstep:respwindowms(2)-cumsumstep ,cumsum(histcounts(indssfirstlicks,respwindowms(1):cumsumstep:respwindowms(2),'Normalization','count'))/numel(indssfirstlicks),'k','LineWidth',2)
            xlim([0 1000])
        end
        if centers(ldelay)<0
            %centers(ldelay)
            hold on;plot(respwindowms(1):cumsumstep:respwindowms(2)-cumsumstep ,cumsum(histcounts(indssfirstlicks,respwindowms(1):cumsumstep:respwindowms(2),'Normalization','count'))/numel(indssfirstlicks),'r','LineWidth',1)
            xlim([0 1000])
        end
        % end
    end
end
indss0= intersect(find(isnan(LaserDelayBinned)),nogotrialind_w);
length(indss0);
crrate(1) = length(intersect(correctnogotrialind_w,indss0))/length(indss0);
farate(1) = length(intersect(incorrectnogotrialind_w,indss0))/length(indss0);
indssfirstlicks=firstlicksample(indss0)/30;

indssfirstlicks(find(indssfirstlicks<respwindowms(1))) =nan;
indssfirstlicks(find(indssfirstlicks>=respwindowms(2))) =nan;
if normalizeforcorrect
    indssfirstlicks(find(isnan(indssfirstlicks))) = [];
end
hold on;plot(respwindowms(1):cumsumstep:respwindowms(2)-cumsumstep ,cumsum(histcounts(indssfirstlicks,respwindowms(1):cumsumstep:respwindowms(2),'Normalization','count'))/numel(indssfirstlicks),'g','LineWidth',1)
xlim([0 1000])

subplot(6,2,6);p3=plot(xpoints,crrate,'+-','Color',[1 0 0]);
p3.Parent.YAxis.Label.String = 'cr rate';
p3.Parent.YAxis.Label.FontWeight='normal';
p3.Parent.YAxis.Label.FontSize=7;
p3.Parent.XLim = [-50 500];
line(p3.Parent.XLim,[crrate(1) crrate(1)],'Color',[0.8 0.8 0.8]);

subplot(6,2,8);p4=plot(xpoints,farate,'+-','Color',[0.5 0.2 0.2]);
p4.Parent.YAxis.Label.String = 'fa rate';
p4.Parent.YAxis.Label.FontWeight='normal';
p4.Parent.YAxis.Label.FontSize=7;
p4.Parent.XLim = [-50 500];
line(p4.Parent.XLim,[farate(1) farate(1)],'Color',[0.8 0.8 0.8]);
%%%%% go and nogo combined
%%%%% all correct trials
correctrate = nan;
incorrectrate = nan;
for ldelay =1:length(centers)    
     indss= find(LaserDelayBinned == (ldelay));
    if length(indss) > min_n_trials_per_delay       
        correctrate(end+1) = length(intersect([correctnogotrialind_w correctgotrialind_w],indss))/length(indss);
        incorrectrate(end+1) = length(intersect([incorrectnogotrialind_w incorrectgotrialind_w],indss))/length(indss);
     end
end
indss0= intersect(find(isnan(LaserDelayBinned)),[gotrialind_w;nogotrialind_w]);
correctrate(1) = length(intersect([correctnogotrialind_w correctgotrialind_w],indss0))/length(indss0);
incorrectrate(1) = length(intersect([incorrectnogotrialind_w incorrectgotrialind_w],indss0))/length(indss0);

subplot(6,2,10);p5=plot(xpoints,correctrate,'+-','Color',[0 0 0]);
p5.Parent.YAxis.Label.String = 'correct rate';
p5.Parent.YAxis.Label.FontWeight='normal';
p5.Parent.YAxis.Label.FontSize=7;
p5.Parent.XLim = [-50 500];
line(p5.Parent.XLim,[correctrate(1) correctrate(1)],'Color',[0.8 0.8 0.8]);

subplot(6,2,12);p6=plot(xpoints,incorrectrate,'+-','Color',[0.8 0.8 0.8]);
p6.Parent.YAxis.Label.String = 'incorrect rate';
p6.Parent.YAxis.Label.FontWeight='normal';
p6.Parent.YAxis.Label.FontSize=7;
p6.Parent.XLim = [-50 500];
line(p6.Parent.XLim,[incorrectrate(1) incorrectrate(1)],'Color',[0.8 0.8 0.8]);