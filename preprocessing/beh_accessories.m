figure;imagesc(licks(incorrectnogotrialind,:))
%%
f_lick2=figure;histogram(firstlicksample(correctgotrialind)/30, 100,'Normalization','count','FaceColor','b')
hold on;histogram(firstvalvesample(correctgotrialind)/30, 100,'Normalization','count','FaceColor','r')
%%
f_lick3=figure;subplot(1,2,1);
histogram(firstvalvesample(correctgotrialind)/30-firstlicksample(correctgotrialind)/30, ...
    100,'Normalization','count','FaceColor','b')
subplot(1,2,2);
histogram(firstvalvesample(intersect(correctgotrialind,nogroomingind))/30-firstlicksample(intersect(correctgotrialind,nogroomingind))/30, ...
    100,'Normalization','count','FaceColor','b')

%%
figure;plot((1:60001)/30 - 1000,mean(licks(intersect(correctgotrialind,nogroomingind),:),1),'Color',[0 1 0])
hold on;plot((1:60001)/30 - 1000,mean(licks(intersect(incorrectgotrialind,nogroomingind),:),1),'Color',[0 0.5 0])
hold on;plot((1:60001)/30 - 1000,mean(licks(intersect(correctnogotrialind,nogroomingind),:),1),'Color',[1 0 0])
hold on;plot((1:60001)/30 - 1000,mean(licks(intersect(incorrectnogotrialind,nogroomingind),:),1),'Color',[0.5 0 0])
%% d-prime

H = length(correctgotrialind)/(length(correctgotrialind)+length(incorrectgotrialind));

F = length(incorrectnogotrialind)/(length(correctnogotrialind)+length(incorrectnogotrialind));


dp=norminv(H) - norminv(F)


H = length(intersect(correctgotrialind,nogroomingind))/...
    (length(intersect(correctgotrialind,nogroomingind))+length(intersect(incorrectgotrialind,nogroomingind)));

F = length(intersect(incorrectnogotrialind,nogroomingind))/...
    (length(intersect(correctnogotrialind,nogroomingind))+length(intersect(incorrectnogotrialind,nogroomingind)));

dp_nogrooming=norminv(H) - norminv(F)