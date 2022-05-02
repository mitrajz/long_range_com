function LAvchdata = fillindedchanels(LAvchdata,dedchs)
for i = 1:length(dedchs)
    i_above = dedchs(i)-1;
    while numel(intersect(dedchs,i_above))
        i_above = i_above -1;
    end
    i_below = dedchs(i)+1;
    while numel(intersect(dedchs,i_below))
        i_below = i_below + 1;
    end
    if dedchs(i) == 1
        LAvchdata(dedchs(i),:) = LAvchdata(i_below,:);
    elseif dedchs(i) == size(LAvchdata,1)
        LAvchdata(dedchs(i),:) = LAvchdata(i_above,:);
    else
        LAvchdata(dedchs(i),:) = (LAvchdata(i_below,:)+LAvchdata(i_above,:))/2;
    end   
end
end

