function [taugo,taunogo,taudif] = expstats(matgo,matnogo,ftype,timescalems)

times = repmat(0:7,size(matgo,1),1);
xgo=reshape(times,1,[])';
ygo=reshape(matgo,1,[])';
novals=isnan(ygo);
xgo(novals)=[];ygo(novals)=[];
xgo = timescalems*xgo;

times = repmat(0:7,size(matnogo,1),1);
xnogo=reshape(times,1,[])';
ynogo=reshape(matnogo,1,[])';
novals=isnan(ynogo);
xnogo(novals)=[];ynogo(novals)=[];
xnogo = timescalems*xnogo;


[mdlgo,~] = fit(xgo,ygo,ftype);
[mdlnogo,~] = fit(xnogo,ynogo,ftype);

taugo.value = mdlgo.tau;
taunogo.value = mdlnogo.tau;

diftau = mdlgo.tau - mdlnogo.tau;
taudif.value = diftau;

%%%% making _lim
fo_lim = fitoptions('Method','NonlinearLeastSquares','Lower',[0],'Upper',[inf],'StartPoint',[1000],'MaxIter',1000);
go_ftype_lim = fittype(sprintf('%d*(exp(-x/tau)+%d)',mdlgo.A,mdlgo.B),'independent','x','options',fo_lim);
nogo_ftype_lim = fittype(sprintf('%d*(exp(-x/tau)+%d)',mdlnogo.A,mdlnogo.B),'independent','x','options',fo_lim);

go_fitdatahandle = @(data)fitdata(data,go_ftype_lim);
nogo_fitdatahandle = @(data)fitdata(data,nogo_ftype_lim);
%fitdifdatahandle = @(datago,datanogo)fitdata_dif(datago,datanogo,ftype);

rng(1,'twister');
[bcigo,statssgo]=bootci(100,go_fitdatahandle,[xgo,ygo]);
rng(1,'twister');
[bcinogo,statssnogo]=bootci(100,nogo_fitdatahandle,[xnogo,ynogo]);

taugo.ci = bcigo;
taunogo.ci = bcinogo;


pooled = [[xgo,ygo];[xnogo,ynogo]];
rng(1,'twister');
tauperm = [];
for i = 1:1000
    ind1 = randsample(size(pooled,1),round(size(pooled,1)/2));
    ind2 = randsample(size(pooled,1),round(size(pooled,1)/2));
    tauperm(end+1) = fitdata_dif(pooled(ind1,:),pooled(ind2,:),ftype);
end
taudif.perm = tauperm;
% two sided test from one sided
numel(find(tauperm>diftau))/numel(tauperm);
numel(find(tauperm<diftau))/numel(tauperm);
taudif.pval = normcdf(-abs(diftau),mean(tauperm),std(tauperm)) + (1-(normcdf(abs(diftau),mean(tauperm),std(tauperm))));

end
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
function out = fitdata_dif(datago,datanogo,ftype)
xgo = datago(:,1);
ygo = datago(:,2);
xnogo = datanogo(:,1);
ynogo = datanogo(:,2);
try
    [mdlgo,~] = fit(xgo,ygo,ftype);
    [mdlnogo,~] = fit(xnogo,ynogo,ftype);
    out = mdlgo.tau - mdlnogo.tau;
catch
    %disp('no fit in jackknife estimate, insufficient data')
    out = nan;
end
end
