function [tau,gof] = cor_timeconstant_plot_helper(mat,color,scAlpha,xjitter,sc,ftype,fo,timescalems,fitmdl)

findbeststart = 0; % default: doesn't optimize starting point, optimization is extremely time taking
                   % and for communication, doesn't really change the fit
                   % result

%%% assign lag1/lag0 values to tau (although not a fit value but raw data)
tau.lag0 = mat(:,1);
tau.lag1 = mat(:,2);
tau.lag8 = mat(:,8);

% prepare for the fit
times = repmat(0:7,size(mat,1),1);
x=reshape(times,1,[])';
y=reshape(mat,1,[])';
novals=isnan(y);
x(novals)=[];y(novals)=[];

% average interval x = k (i-j)
x = timescalems*x;

if strcmp(fitmdl,'exp')
    % find the best start point: the point that minimizes square error is
    % chosen and used as fo initial point. The same fo is then chosen in the
    % jackknife. Time taking process, comment out for a quick fit
    try
        if findbeststart
            [bestmdl,beststart] = findbeststart_forfit(x,y,fo,ftype);
            fo.StartPoint= beststart;
            ftype = setoptions(ftype,fo);
        else
            bestmdl.tau = nan;
        end
        
        
        % model fitted to the full data. tau estimate calculated
        [mdl,gof] = fit(x,y,ftype);
        %mdl = bestmdl;gof=nan;
        [mdl.tau bestmdl.tau]
        tau.val = mdl.tau;
        tau.bestmdlvalue = bestmdl.tau;
        
        % jackknife estimate of standard error (95% ci), bias not corrected (large
        % n)
        fitdatahandle = @(data)fitdata(data,ftype);
        tau_j = jackknife(fitdatahandle,[x,y]);
        % nanmean(tau_j) This should be very close to tau.val (checked)
        tau.ci = 2*nanstd(tau_j);
        
        %%% plot:
        
        hold on;
        plot(unique(x),mdl.A*(exp(-unique(x)/mdl.tau)+mdl.B),'-','Color',color)
        hold on;
        
        
        yp = predint(mdl,unique(x),0.95,'functional','on');
        patch([unique(x)',fliplr(unique(x)')],...
            [yp(:,1)',fliplr(yp(:,2)')],color,...
            'EdgeAlpha',0,'FaceAlpha',0.2);
        
        hold on;
        if sc
            rndwdth = 0.1;
            xp = x - rndwdth/2 + rndwdth*rand(size(x));
            scatter((xp+xjitter/timescalems),y,'.','MarkerEdgeColor',color,'MarkerEdgeAlpha',scAlpha)
            hold on;
        end
        ylim([0,1])
        
        
        if strcmp(color,'r')|strcmp(color,'k')
            text(1.5,0.05,sprintf('tau = %f(ci:%f), A = %f, B = %f',mdl.tau,tau.ci,mdl.A,mdl.B),'Color',color)
        elseif strcmp(color,'g')|strcmp(color,'c')
            text(1.5 ,0.1,sprintf('tau = %f(ci:%f), A = %f, B = %f',mdl.tau,tau.ci,mdl.A,mdl.B),'Color',color)
        end
    catch
        disp('not enough data for tc plots/confidence interval')
        tau.val = nan;
        tau.ci = nan;
        gof = nan;
    end
elseif strcmp(fitmdl , 'lin')
    [mdl,gof] = fit(x,y,ftype);
    % in ta slope (A) is stored instead of tau in exp case
    tau.val = mdl.A;
    tau.ci = nan;
    
    fitdatahandle_lin = @(data)fitdata_lin(data,ftype);
    tau_j = jackknife(fitdatahandle_lin,[x,y]);
    tau.ci = 2*nanstd(tau_j);
    
    % plots
    hold on;
    plot(unique(x),mdl.A*unique(x)+mdl.B,'-','Color',color)
    hold on;
    
     
    yp = predint(mdl,unique(x),0.95,'functional','on');
    patch([unique(x)',fliplr(unique(x)')],...
        [yp(:,1)',fliplr(yp(:,2)')],color,...
        'EdgeAlpha',0,'FaceAlpha',0.2);
    
    hold on;
    if sc
        rndwdth = 0.1;
        xp = x - rndwdth/2 + rndwdth*rand(size(x));
        scatter((xp+xjitter/timescalems),y,'.','MarkerEdgeColor',color,'MarkerEdgeAlpha',scAlpha)
        hold on;
    end
    ylim([0,1])
    
    if strcmp(color,'r')|strcmp(color,'k')
        text(1.5,0.05,sprintf('A = %f(ci:%f), B = %f',mdl.A,tau.ci,mdl.B),'Color',color)
    elseif strcmp(color,'g')|strcmp(color,'c')
        text(1.5 ,0.1,sprintf('A = %f(ci:%f), B = %f',mdl.A,tau.ci,mdl.B),'Color',color)
    end
else
    disp('skipped expfit')
    tau.val = nan;
    tau.ci = nan;
    gof = nan;
    
end
% plot mean and 95% ci of mean (2 ste) for each x, over all points
mat_mean = nanmean(mat,1);
mat_std = nanstd(mat);
mat_num_elements = sum(~isnan(mat),1);
mat_ste = mat_std./sqrt(mat_num_elements);

mat_mean = mat_mean(~isnan(mat_mean));
mat_ste = mat_ste(~isnan(mat_ste));

errorbar(xjitter+unique(x),mat_mean,2*mat_ste,'.','Color',color,'LineStyle','none')


end
%%
function out = fitdata(data,ftype)
x = data(:,1);
y = data(:,2);
try
    [mdl,gof] = fit(x,y,ftype);
    out = mdl.tau;
catch
    %disp('no fit in jackknife estimate, insufficient data')
    out = nan;
end
end
function out = fitdata_lin(data,ftype)
x = data(:,1);
y = data(:,2);
try
    [mdl,gof] = fit(x,y,ftype);
    out = mdl.A;
catch
    %disp('no fit in jackknife estimate, insufficient data')
    out = nan;
end
end