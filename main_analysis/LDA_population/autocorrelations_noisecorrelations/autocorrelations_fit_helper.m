function [mdl,gof] = autocorrelations_fit_helper(cormat_tc_go,ftype,timescalems)

mat = cormat_tc_go(:,2:end);
times = repmat(1:7,size(mat,1),1);
x=reshape(times,1,[])';
y=reshape(mat,1,[])';
novals=isnan(y)
x(novals)=[];y(novals)=[];
if numel(x)<5
    mdl = nan;gof = nan;
    
else
    [mdl,gof] = fit(x,y,ftype);
    
%     figure;
%     hold on;
%     plot(timescalems*(1:7),mdl.A*(exp(mdl.b*(1:7))+mdl.D),'-','Color','g')
%     hold on;
%     yp = predint(mdl,1:7,0.99,'functional','on');
%     patch(timescalems*[1:7,7:-1:1],...
%         [yp(:,1)',fliplr(yp(:,2)')],'g',...
%         'EdgeAlpha',0,'FaceAlpha',0.2);
    
    hold on;
    
   % scatter(timescalems*(x),y,'.','MarkerEdgeColor','g','MarkerEdgeAlpha',1)
    hold on;
    title(sprintf('%d',gof.rmse))
end