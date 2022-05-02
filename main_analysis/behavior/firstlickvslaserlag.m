
figure;
normalizeforcorrect = 0;
for ldelay =1:length(centers)  
  indss= intersect(find(LaserDelayBinned == (ldelay)),[gotrialind]);     


    if length(indss) > min_n_trials_per_delay
        length(indss)
        indssfirstlicks=firstlicksample(indss)/30;
     
       indssfirstlicks(find(indssfirstlicks<150)) =nan;
       indssfirstlicks(find(indssfirstlicks>=690)) =nan;
       if normalizeforcorrect
            indssfirstlicks(find(isnan(indssfirstlicks))) = [];
       end
        hold on; plot(100:20:680-1,cumsum(histcounts(indssfirstlicks,100:20:680,'Normalization','count'))/numel(indssfirstlicks),'b')
        xlim([0 1000])

        if centers(ldelay)> 50 & centers(ldelay)< 60
            centers(ldelay)
            hold on;plot(100:20:680-1 ,cumsum(histcounts(indssfirstlicks,100:20:680,'Normalization','count'))/numel(indssfirstlicks),'k','LineWidth',2)
            xlim([0 1000])
        end
        if centers(ldelay)<0
            centers(ldelay)
            hold on;plot(100:20:680-1 ,cumsum(histcounts(indssfirstlicks,100:20:680,'Normalization','count'))/numel(indssfirstlicks),'r','LineWidth',1)
            xlim([0 1000])
        end

    end
    
    indss0= intersect(find(isnan(LaserDelayBinned)),[gotrialind]);    
    length(indss0)
    indssfirstlicks=firstlicksample(indss0)/30;
    
    indssfirstlicks(find(indssfirstlicks<150)) =nan;
    indssfirstlicks(find(indssfirstlicks>=690)) =nan;
    if normalizeforcorrect
        indssfirstlicks(find(isnan(indssfirstlicks))) = [];
    end
    hold on;plot(100:20:680-1 ,cumsum(histcounts(indssfirstlicks,100:20:680,'Normalization','count'))/numel(indssfirstlicks),'g','LineWidth',1)
    xlim([0 1000])
end

