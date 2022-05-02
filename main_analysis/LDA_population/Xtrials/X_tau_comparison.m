
%% params
ldaormu = 1; % if 1, all for lda if 2: mu but  when 2 importantly
preSub = ''; % options: '' (default), '_preSub', and '_preSub2'
cidivide = 1; % if 2:sem, if 1:95%ci
% activity from 65 and com from 150
% these are results from X_cor_timeconstant.m
root = '/mnt/data/Mitra/figs/P2_L/bothdircombined/2020_02_24/LDA/version2/cosinesimilaritymats_and_timeconstants/saveres/';
% tau: com vs source activity per animal
cd(root);
if strcmp(preSub,'')
    cd('FF150Xstyle2nrep100_withsource_2020_11_09_21_05_29');
    %cd('FF150Xstyle2nrep100_withsource_2020_11_05_15_22_55');
elseif strcmp(preSub,'_preSub')
    cd('FF150Xstyle2nrep100_withsource_ldaormu1_preSub_2021_12_06_13_00_35');
elseif strcmp(preSub,'_preSub2')
    cd('FF150Xstyle2nrep100_withsource_ldaormu1_preSub2_2021_12_06_16_15_55');
end
ff150 = load('alltaus.mat');


cd(root);
if strcmp(preSub,'')
    cd('notau/FF65Xstyle2nrep100_withsource_2020_11_09_20_29_35');
    %cd('FF65Xstyle2nrep100_withsource_2020_11_05_15_15_48');
elseif strcmp(preSub,'_preSub')
    cd('FF65Xstyle2nrep100_withsource_ldaormu1_preSub_2021_12_06_12_54_10');
elseif strcmp(preSub,'_preSub2')
    cd('FF65Xstyle2nrep100_withsource_ldaormu1_preSub2_2021_12_06_16_09_34');
end
ff65 = load('alltaus.mat');

cd(root);
if strcmp(preSub,'')
    cd('FB150Xstyle2nrep100_withsource_2020_11_09_21_10_39');
    %cd('FB150Xstyle2nrep100_withsource_2020_11_05_15_35_26');
elseif strcmp(preSub,'_preSub')
    cd('FB150Xstyle2nrep100_withsource_ldaormu1_preSub_2021_12_06_13_05_04');
elseif strcmp(preSub,'_preSub2')
    cd('FB150Xstyle2nrep100_withsource_ldaormu1_preSub2_2021_12_06_16_23_48');
end
fb150 = load('alltaus.mat');

cd(root);
if strcmp(preSub,'')
    cd('notau/FB65Xstyle2nrep100_withsource_2020_11_09_20_30_35');
    %cd('FB65Xstyle2nrep100_withsource_2020_11_05_15_28_48');
elseif strcmp(preSub,'_preSub')
    cd('FB65Xstyle2nrep100_withsource_ldaormu1_preSub_2021_12_06_12_46_53');
elseif strcmp(preSub,'_preSub2')
    cd('FB65Xstyle2nrep100_withsource_ldaormu1_preSub2_2021_12_06_16_19_17');
end
fb65 = load('alltaus.mat');

cd(root);
cd('notau/simul_FB150_2020_11_09_20_32_33');
%cd('FB65Xstyle2nrep100_withsource_2020_11_05_15_28_48');
fbs = load('alltaus.mat');

cd(root);
cd('notau/simul_FF150_2020_11_09_20_31_55');
%cd('FB65Xstyle2nrep100_withsource_2020_11_05_15_28_48');
ffs = load('alltaus.mat');

% in untrained animals go = nogo = average of the 2 stimuli
cd(root);
cd('Untrained-FB150Xstyle2nrep100_withsource_2020_11_16_16_37_44');
%cd('FB65Xstyle2nrep100_withsource_2020_11_05_15_28_48');
fbu = load('alltaus.mat');

% pc slopes
cd(root);
cd('pc_FB150Xstyle1nrep100_eqbsls4numtrials10_2020_12_29_00_50_18');
pc = load('allanimalpcs.mat');

% fast and slow cells
cd(root);
cd('V1-slow0.8_Activity_ff_fb_combined65Xstyle2nrep100_2021_01_10_22_43_31');
%cd('V1-slow0.8_Activity_ff_fb_combined65Xstyle2nrep100_2020_11_18_16_39_10');
slowV1 = load('alltaus.mat');

cd(root);
cd('V1-fast0.2_Activity_ff_fb_combined65Xstyle2nrep100_2021_01_10_22_45_43');
%cd('V1-fast0.2_Activity_ff_fb_combined65Xstyle2nrep100_2020_11_18_16_35_56');
fastV1 = load('alltaus.mat');

cd(root);
cd('LM-slow0.8_Activity_ff_fb_combined65Xstyle2nrep100_2021_01_10_22_41_32');
%cd('LM-slow0.8_Activity_ff_fb_combined65Xstyle2nrep100_2020_11_18_16_37_29');
slowLM = load('alltaus.mat');

cd(root);
cd('LM-fast0.2_Activity_ff_fb_combined65Xstyle2nrep100_2021_01_10_22_39_16');
%cd('LM-fast0.2_Activity_ff_fb_combined65Xstyle2nrep100_2020_11_18_16_34_15');
fastLM = load('alltaus.mat');

cd(root)
cd('FF150Xstyle2nrep100_withsource_ldaormu2_2021_01_04_18_06_31');
ffmu150 = load('alltaus.mat');

cd(root)
cd('FB150Xstyle2nrep100_withsource_ldaormu2_2021_01_04_19_34_34')
fbmu150 = load('alltaus.mat');

cd(root);
cd('targetreduced20FF150Xstyle2nrep100_withsource_ldaormu1_2021_02_20_00_06_04');
ffreduced = load('alltaus.mat');

cd(root);
cd('targetreduced20FB150Xstyle2nrep100_withsource_ldaormu1_2021_02_20_00_08_11');
fbreduced = load('alltaus.mat');


if ldaormu == 2
    fb150 = fbmu150;
    ff150 = ffmu150;
end


%%% ff
gotarget = cellfun(@(x)x.val, ff65.animaltau_n0_go);
nogotarget = cellfun(@(x)x.val, ff65.animaltau_n0_nogo);
gocom = cellfun(@(x)x.val, ff150.animaltau_com_go);
nogocom = cellfun(@(x)x.val, ff150.animaltau_com_nogo);

ci_gotarget = cellfun(@(x)x.ci, ff65.animaltau_n0_go)/cidivide;
ci_nogotarget = cellfun(@(x)x.ci, ff65.animaltau_n0_nogo)/cidivide;
ci_gocom = cellfun(@(x)x.ci, ff150.animaltau_com_go)/cidivide;
ci_nogocom = cellfun(@(x)x.ci, ff150.animaltau_com_nogo)/cidivide;

figure;plotgonogopairs(gotarget,gocom,nogotarget,nogocom,...
    ci_gotarget,ci_gocom,ci_nogotarget,ci_nogocom)


%%% fb

gotarget = cellfun(@(x)x.val, fb65.animaltau_n0_go);
nogotarget = cellfun(@(x)x.val, fb65.animaltau_n0_nogo);
gocom = cellfun(@(x)x.val, fb150.animaltau_com_go);
nogocom = cellfun(@(x)x.val, fb150.animaltau_com_nogo);

ci_gotarget = cellfun(@(x)x.ci, fb65.animaltau_n0_go)/cidivide;
ci_nogotarget = cellfun(@(x)x.ci, fb65.animaltau_n0_nogo)/cidivide;
ci_gocom = cellfun(@(x)x.ci, fb150.animaltau_com_go)/cidivide;
ci_nogocom = cellfun(@(x)x.ci, fb150.animaltau_com_nogo)/cidivide;

figure;plotgonogopairs(gotarget,gocom,nogotarget,nogocom,...
    ci_gotarget,ci_gocom,ci_nogotarget,ci_nogocom)

%%%%% allanimalstogether
gotarget = cellfun(@(x)x.val, ff65.tau_n0_go);
nogotarget = cellfun(@(x)x.val, ff65.tau_n0_nogo);
gocom = cellfun(@(x)x.val, ff150.tau_com_go);
nogocom = cellfun(@(x)x.val, ff150.tau_com_nogo);

ci_gotarget = cellfun(@(x)x.ci, ff65.tau_n0_go)/cidivide;
ci_nogotarget = cellfun(@(x)x.ci, ff65.tau_n0_nogo)/cidivide;
ci_gocom = cellfun(@(x)x.ci, ff150.tau_com_go)/cidivide;
ci_nogocom = cellfun(@(x)x.ci, ff150.tau_com_nogo)/cidivide;

figure;plotgonogopairs(gotarget,gocom,nogotarget,nogocom,...
    ci_gotarget,ci_gocom,ci_nogotarget,ci_nogocom)

gotarget = cellfun(@(x)x.val, fb65.tau_n0_go);
nogotarget = cellfun(@(x)x.val, fb65.tau_n0_nogo);
gocom = cellfun(@(x)x.val, fb150.tau_com_go);
nogocom = cellfun(@(x)x.val, fb150.tau_com_nogo);

ci_gotarget = cellfun(@(x)x.ci, fb65.tau_n0_go)/cidivide;
ci_nogotarget = cellfun(@(x)x.ci, fb65.tau_n0_nogo)/cidivide;
ci_gocom = cellfun(@(x)x.ci, fb150.tau_com_go)/cidivide;
ci_nogocom = cellfun(@(x)x.ci, fb150.tau_com_nogo)/cidivide;

hold on;plotgonogopairs(gotarget,gocom,nogotarget,nogocom,...
    ci_gotarget,ci_gocom,ci_nogotarget,ci_nogocom)
%% tau single animal - go com an nogo com in go and nogo

yl = [0 1000];
figure;
subplot(1,2,1);
hold on;multiplelines(cellfun(@(x) x.val, ff150.animaltau_com_go),...
    cellfun(@(x) x.val, ff150.animaltau_com_nogo),1,2,yl)
hold on;multiplelines(cellfun(@(x) x.val, fb150.animaltau_com_go),...
    cellfun(@(x) x.val, fb150.animaltau_com_nogo),3,4,yl)

xticks([1 2 3 4])
xticklabels({'ffcom-go','ffcom-nogo','fbcom-go','fbcom-nogo'})
xlim([0 5]); ylim([yl(1),yl(2)+200])
ylabel('tau')

subplot(1,2,2);
hold on;multiplebars(cellfun(@(x) x.val, ff150.animaltau_com_go),...
    cellfun(@(x) x.val, ff150.animaltau_com_nogo),1,2,yl,'g','r')
hold on;multiplebars(cellfun(@(x) x.val, fb150.animaltau_com_go),...
    cellfun(@(x) x.val, fb150.animaltau_com_nogo),3,4,yl,'g','r')
% tau estimate doesn't work for untrained,  3 animals, lots of variability
hold on;u_multiplebars(cellfun(@(x) x.val, fbu.animaltau_com_go),...
    cellfun(@(x) x.val, fb150.animaltau_com_go),...
    cellfun(@(x) x.val, fb150.animaltau_com_nogo),5,yl,'k')

xticks([1 2 3 4])
xticklabels({'ffcom-go','ffcom-nogo','fbcom-go','fbcom-nogo'})
xlim([0 6]); ylim([yl(1),yl(2)+200])
ylabel('tau')
%% lago,lag1,lag 8 (slope and offset)
% in scatter plots ffgo, fbgo, ffnogo, and fbnogo are scattered separately,
% so can be plotted by 4 different colors, or markers. However, for the
% current version, go and nogo might be combined
fastcom=figure;


op = 'slope'; % slope, norml1, norml8
type = 'peranimal';% 'peranimal', or 'perlag' : per animal for most things(Q1), per lag for slope go/nogo comparison (Q2)
deltat = 0.065; % perlag are usually the bar types, per animal are scatters 
propname = op;

if strcmp(op,'slope')
yl = [-10,0];
else
yl = [-0.5,1.5];
end
%%%%% lag1/lag0
% this types: ranksum or signrank
% com slower than conneciviy in both ff and fb
% 1- com vs target
figure(fastcom);subplot(1,3,1);line(yl,yl,'Color','k');
% ff target vs ff com go
hold on;scatter(getallvals(ff65.animaltau_n0_go,op,type,deltat),...
    getallvals(ff150.animaltau_com_go,op,type,deltat),200,[1 0.65 0],'.')
% ff target vs ff com nogo
hold on;scatter(getallvals(ff65.animaltau_n0_nogo,op,type,deltat),...
    getallvals(ff150.animaltau_com_nogo,op,type,deltat),200,[1 0.65 0],'.')
% fb target vs fb com go
hold on;scatter(getallvals(fb65.animaltau_n0_go,op,type,deltat),...
    getallvals(fb150.animaltau_com_go,op,type,deltat),200,[0.6 0.11 0.72],'.')
% fb target vs fb com nogo
hold on;scatter(getallvals(fb65.animaltau_n0_nogo,op,type,deltat),...
    getallvals(fb150.animaltau_com_nogo,op,type,deltat),200,[0.6 0.11 0.72],'.')
xlabel(['target',propname]); ylabel(['com',propname])
legend({'','ff-go','ff-nogo','fb-go','fb-nogo'},'location','northwest')

p_ov=signrank([getallvals(ff65.animaltau_n0_go,op,type,deltat),getallvals(ff65.animaltau_n0_nogo,op,type,deltat),...
    getallvals(fb65.animaltau_n0_go,op,type,deltat),getallvals(fb65.animaltau_n0_nogo,op,type,deltat)],...
    [getallvals(ff150.animaltau_com_go,op,type,deltat),getallvals(ff150.animaltau_com_nogo,op,type,deltat),...
    getallvals(fb150.animaltau_com_go,op,type,deltat),getallvals(fb150.animaltau_com_nogo,op,type,deltat)]);
hold on; text(-6,0,['pval signrank = ',num2str(p_ov)])

% 2- com vs source: make sure n1 is source no normalized target
figure(fastcom);subplot(1,3,2);line(yl,yl,'Color','k');
% ff source vs ff com go
hold on;scatter(getallvals(ff65.animaltau_n1_go,op,type,deltat),...
    getallvals(ff150.animaltau_com_go,op,type,deltat),200,[1 0.65 0],'.')
% ff source vs ff com nogo
hold on;scatter(getallvals(ff65.animaltau_n1_nogo,op,type,deltat),...
    getallvals(ff150.animaltau_com_nogo,op,type,deltat),200,[1 0.65 0],'.')
% fb source vs fb com go
hold on;scatter(getallvals(fb65.animaltau_n1_go,op,type,deltat),...
    getallvals(fb150.animaltau_com_go,op,type,deltat),200,[0.6 0.11 0.72],'.')
% fb source vs fb com nogo
hold on;scatter(getallvals(fb65.animaltau_n1_nogo,op,type,deltat),...
    getallvals(fb150.animaltau_com_nogo,op,type,deltat),200,[0.6 0.11 0.72],'.')

xlabel(['source',propname]); ylabel(['com',propname])
legend({'','ff-go','ff-nogo','fb-go','fb-nogo'},'location','northwest')

p_ov=signrank([getallvals(ff65.animaltau_n1_go,op,type,deltat),getallvals(ff65.animaltau_n1_nogo,op,type,deltat),...
    getallvals(fb65.animaltau_n1_go,op,type,deltat),getallvals(fb65.animaltau_n1_nogo,op,type,deltat)],...
    [getallvals(ff150.animaltau_com_go,op,type,deltat),getallvals(ff150.animaltau_com_nogo,op,type,deltat),...
    getallvals(fb150.animaltau_com_go,op,type,deltat),getallvals(fb150.animaltau_com_nogo,op,type,deltat)]);
hold on; text(-6,0,['pval signrank = ',num2str(p_ov)])


% 3- com: ff slope go = slope nogo, fb slope go = slope nogo
figure;line(yl,yl,'Color','k');
% ff source vs ff com go
hold on;scatter(getallvals(ff150.animaltau_com_nogo,op,type,deltat),...
    getallvals(ff150.animaltau_com_go,op,type,deltat),100,'k.')
xlabel(['nogo com',propname]);ylabel(['go com',propname]);
% fb source vs fb com go
hold on;scatter(getallvals(fb150.animaltau_com_nogo,op,type,deltat),...
    getallvals(fb150.animaltau_com_go,op,type,deltat),100,'k*')
xlabel(['nogo com',propname]);ylabel(['go com',propname]);
legend({'','ff','fb'})

% 4 - same as above but in 2 columns
figure;
subplot(1,2,1);
hold on;multiplelines(getallvals(ff150.animaltau_com_go,op,type,deltat),...
    getallvals(ff150.animaltau_com_nogo,op,type,deltat),1,2,yl)
hold on;multiplelines(getallvals(fb150.animaltau_com_go,op,type,deltat),...
    getallvals(fb150.animaltau_com_nogo,op,type,deltat),3,4,yl)

xticks([1 2 3 4])
xticklabels({'ffcom-go','ffcom-nogo','fbcom-go','fbcom-nogo'})

ylabel(propname)
xlim([0 5])
ylim([yl(1),yl(2)+0.6])

subplot(1,2,2);
hold on;multiplebars(getallvals(ff150.animaltau_com_go,op,type,deltat),...
    getallvals(ff150.animaltau_com_nogo,op,type,deltat),1,2,yl,'g','r')
hold on;multiplebars(getallvals(fb150.animaltau_com_go,op,type,deltat),...
    getallvals(fb150.animaltau_com_nogo,op,type,deltat),3,4,yl,'g','r')

hold on;u_multiplebars(getallvals(fbu.animaltau_com_go,op,type,deltat),...
getallvals(fb150.animaltau_com_go,op,type,deltat),...
getallvals(fb150.animaltau_com_nogo,op,type,deltat),5,yl,'k')


xticks([1 2 3 4])
xticklabels({'ffcom-go','ffcom-nogo','fbcom-go','fbcom-nogo'})

ylabel(propname)
xlim([0 6])
ylim([yl(1),yl(2)+0.6])

%  5 - com vs simul com
% ff simul vs ff com go
figure(fastcom);subplot(1,3,3);line(yl,yl,'Color','k');
hold on;scatter(getallvals(ffs.animaltau_com_go,op,type,deltat),...
    getallvals(ff150.animaltau_com_go,op,type,deltat),200,[1 0.65 0],'.')
% ff simul vs ff com nogo
hold on;scatter(getallvals(ffs.animaltau_com_nogo,op,type,deltat),...
    getallvals(ff150.animaltau_com_nogo,op,type,deltat),200,[1 0.65 0],'.')
% fb simul vs fb com go
hold on;scatter(getallvals(fbs.animaltau_com_go,op,type,deltat),...
    getallvals(fb150.animaltau_com_go,op,type,deltat),200,[0.6 0.11 0.72],'.')
% fb simul vs fb com nogo
hold on;scatter(getallvals(fbs.animaltau_com_nogo,op,type,deltat),...
    getallvals(fb150.animaltau_com_nogo,op,type,deltat),200,[0.6 0.11 0.72],'.')
xlabel(['simul',propname]); ylabel(['com',propname])
legend({'','ff-go','ff-nogo','fb-go','fb-nogo'},'location','northwest')

p_ov=signrank([getallvals(ffs.animaltau_com_go,op,type,deltat),getallvals(ffs.animaltau_com_nogo,op,type,deltat),...
    getallvals(fbs.animaltau_com_go,op,type,deltat),getallvals(fbs.animaltau_com_nogo,op,type,deltat)],...
    [ getallvals(ff150.animaltau_com_go,op,type,deltat),getallvals(ff150.animaltau_com_nogo,op,type,deltat),...
    getallvals(fb150.animaltau_com_go,op,type,deltat),getallvals(fb150.animaltau_com_nogo,op,type,deltat)]);
hold on; text(-6,0,['pval signrank = ',num2str(p_ov)])

% 6- fb: source and target not different in go vs nogo?
figure;
subplot(1,2,1)
hold on;multiplelines(getallvals(fb65.animaltau_n1_go,op,type,deltat),...
    getallvals(fb65.animaltau_n1_nogo,op,type,deltat),1,2,yl)
hold on;multiplelines(getallvals(fb65.animaltau_n0_go,op,type,deltat),...
    getallvals(fb65.animaltau_n0_nogo,op,type,deltat),3,4,yl)

xticks([1 2 3 4])
xticklabels({'fbsourcego','fbsource-nogo','fbtarget-go','fbtarget-nogo'})
ylabel(propname)

xlim([0 5])
ylim([yl(1),yl(2)+0.6])

subplot(1,2,2)
hold on;multiplebars(getallvals(fb65.animaltau_n1_go,op,type,deltat),...
    getallvals(fb65.animaltau_n1_nogo,op,type,deltat),1,2,yl,'g','r')
hold on;multiplebars(getallvals(fb65.animaltau_n0_go,op,type,deltat),...
    getallvals(fb65.animaltau_n0_nogo,op,type,deltat),3,4,yl,'g','r')

xticks([1 2 3 4])
xticklabels({'fbsourcego','fbsource-nogo','fbtarget-go','fbtarget-nogo'})
ylabel(propname)

xlim([0 5])
ylim([yl(1),yl(2)+0.6])

% 7-V1 and LM (ff and fb combined)
figure;
subplot(1,2,2)

hold on;multiplebars([getallvals(fb65.animaltau_n1_go,op,type,deltat),getallvals(ff65.animaltau_n0_go,op,type,deltat)],...
    [getallvals(fb65.animaltau_n1_nogo,op,type,deltat),getallvals(ff65.animaltau_n0_nogo,op,type,deltat)],1,2,yl,'g','r')
hold on;multiplebars([getallvals(fb65.animaltau_n0_go,op,type,deltat),getallvals(ff65.animaltau_n1_go,op,type,deltat)],...
    [getallvals(fb65.animaltau_n0_nogo,op,type,deltat),getallvals(ff65.animaltau_n1_nogo,op,type,deltat)],3,4,yl,'g','r')

xticks([1 2 3 4])
xticklabels({'LM-go','LM-nogo','V1-go','V1-nogo'})
ylabel(propname)

xlim([0 5])
ylim([yl(1),yl(2)+0.6])

%%% example statistic
% ranksum(cellfun(@(x) nanmean(getfield(x,'lag1'))/nanmean(x.lag0), ff65.animaltau_n0_go),...
%     cellfun(@(x) nanmean(getfield(x,'lag1'))/nanmean(x.lag0), ff65.animaltau_n0_nogo))
%

%% PCs
% type = 'perlag' is needed for here
type = 'perlag';
animal = pc.animal;
figure;
yl = [-3,1];
for i =1:3 % differemt pcs
    s=subplot(2,3,i);
    s.Title.String = type;
    st_bs_go = cellfun(@(x) x.pcn{i}.bs.go,animal,'UniformOutput',0);
    st_bs_nogo = cellfun(@(x) x.pcn{i}.bs.nogo,animal,'UniformOutput',0);
    st_ls_go = cellfun(@(x) x.pcn{i}.ls.go,animal,'UniformOutput',0);
    st_ls_nogo = cellfun(@(x) x.pcn{i}.ls.nogo,animal,'UniformOutput',0);
    
    subplot(2,3,i);multiplebars(getallvals(st_bs_go,op,type,deltat),...
        getallvals(st_bs_nogo,op,type,deltat),1,2,yl,'g','r')
    xticks([1 2])
    xticklabels({'bs-go','bs-nogo'})
    ylabel(propname)
    xlim([0 3])
    ylim([yl(1),yl(2)+0.6])
    
    subplot(2,3,3+i);multiplebars(getallvals(st_ls_go,op,type,deltat),...
        getallvals(st_ls_nogo,op,type,deltat),1,2,yl,'g','r')
    xticks([1 2])
    xticklabels({'ls-go','ls-nogo'})
    ylabel(propname)
    xlim([0 3])
    ylim([yl(1),yl(2)+0.6])       
end

st_bs_go = [];
st_bs_nogo = [];
st_ls_go = [];
st_ls_nogo = [];

for i =1:3 % differemt pcs   
    st_bs_go = [st_bs_go,cellfun(@(x) x.pcn{i}.bs.go,animal,'UniformOutput',0)];
    st_bs_nogo = [st_bs_nogo,cellfun(@(x) x.pcn{i}.bs.nogo,animal,'UniformOutput',0)];
    st_ls_go = [st_ls_go,cellfun(@(x) x.pcn{i}.ls.go,animal,'UniformOutput',0)];
    st_ls_nogo = [st_ls_nogo,cellfun(@(x) x.pcn{i}.ls.nogo,animal,'UniformOutput',0)];
end


figure;multiplebars((getallvals(st_bs_go,op,type,deltat)-getallvals(st_bs_nogo,op,type,deltat)),...
    (getallvals(st_ls_go,op,type,deltat)-getallvals(st_ls_nogo,op,type,deltat)),1,2,yl,'k','b')
xticks([1 2])
xticklabels({'baseline','laser'})
ylabel('go-nogo differences, pooled across pcs')
xlim([0 3])
ylim([yl(1),yl(2)+0.6])
%% slow and fast


yl = [-9,2];

%%%%%% V1
type = 'perlag';
figure;
s = subplot(1,3,1);
s.Title.String = type;
hold on;multiplebars(...
    (getallvals(fastV1.animaltau_n0_go,op,type,deltat)+getallvals(fastV1.animaltau_n0_nogo,op,type,deltat))/2,...
    (getallvals(slowV1.animaltau_n0_go,op,type,deltat)+getallvals(slowV1.animaltau_n0_nogo,op,type,deltat))/2,...
    1,2,yl,'m','b')

xticks([1 2])
xticklabels({'fast-v1','slow-v1'})

ylabel(propname)
xlim([0 5])
ylim([yl(1),yl(2)+0.6])

s = subplot(1,3,2);
s.Title.String = type;
hold on;multiplebars(getallvals(fastV1.animaltau_n0_go,op,type,deltat),...
    getallvals(fastV1.animaltau_n0_nogo,op,type,deltat),1,2,yl,'g','r')
hold on;multiplebars(getallvals(slowV1.animaltau_n0_go,op,type,deltat),...
    getallvals(slowV1.animaltau_n0_nogo,op,type,deltat),3,4,yl,'g','r')
xticks([1 2 3 4])
xticklabels({'fast-go-v1','fast-nogo-v1','slow-go-v1','slow-nogo-v1'})

ylabel(propname)
xlim([0 5])
ylim([yl(1),yl(2)+0.6])


type = 'peranimal';
s = subplot(1,3,3);
s.Title.String = type;

line(yl,yl,'Color','k');
% fb fast source (V1) vs fb com go
% order is ff then fb
allanimalsV1go= getallvals(fastV1.animaltau_n0_go,op,type,deltat);
hold on;scatter(allanimalsV1go(1:7),...
    getallvals(ffreduced.animaltau_com_go,op,type,deltat),200,'g','.')
% fb fast source vs fb com nogo
allanimalsV1nogo= getallvals(fastV1.animaltau_n0_nogo,op,type,deltat);
hold on;scatter(allanimalsV1nogo(1:7),...
    getallvals(ffreduced.animaltau_com_nogo,op,type,deltat),200,'r','.')
xlabel('fast V1 activity in ff animals - slope');
ylabel('ff-com-slope');


p_ov=signrank([allanimalsV1go(1:7),allanimalsV1nogo(1:7)],...
    [getallvals(ffreduced.animaltau_com_go,op,type,deltat),getallvals(ffreduced.animaltau_com_nogo,op,type,deltat)]);
hold on; text(-6,0,['pval signrank = ',num2str(p_ov)])

%%%%%% LM
type = 'perlag';
figure;
s=subplot(1,3,1);
s.Title.String = type;
hold on;multiplebars(...
    (getallvals(fastLM.animaltau_n0_go,op,type,deltat)+getallvals(fastLM.animaltau_n0_nogo,op,type,deltat))/2,...
    (getallvals(slowLM.animaltau_n0_go,op,type,deltat)+getallvals(slowLM.animaltau_n0_nogo,op,type,deltat))/2,...
    1,2,yl,'m','b')

xticks([1 2])
xticklabels({'fast-LM','slow-LM'})

ylabel(propname)
xlim([0 5])
ylim([yl(1),yl(2)+0.6])

s = subplot(1,3,2);
s.Title.String = type;
hold on;multiplebars(getallvals(fastLM.animaltau_n0_go,op,type,deltat),...
    getallvals(fastLM.animaltau_n0_nogo,op,type,deltat),1,2,yl,'g','r')
hold on;multiplebars(getallvals(slowLM.animaltau_n0_go,op,type,deltat),...
    getallvals(slowLM.animaltau_n0_nogo,op,type,deltat),3,4,yl,'g','r')
xticks([1 2 3 4])
xticklabels({'fast-go-LM','fast-nogo-LM','slow-go-LM','slow-nogo-LM'})

ylabel(propname)
xlim([0 5])
ylim([yl(1),yl(2)+0.6])

type = 'peranimal';
s = subplot(1,3,3);
s.Title.String = type;
line(yl,yl,'Color','k');
% ff fast source (LM) vs ff com go
% order is fb then ff
allanimalsLMgo= getallvals(fastLM.animaltau_n0_go,op,type,deltat);
hold on;scatter(allanimalsLMgo(1:6),...
    getallvals(fbreduced.animaltau_com_go,op,type,deltat),200,'g','.')
% fb fast source vs fb com nogo
allanimalsLMnogo= getallvals(fastLM.animaltau_n0_nogo,op,type,deltat);
hold on;scatter(allanimalsLMnogo(1:6),...
    getallvals(fbreduced.animaltau_com_nogo,op,type,deltat),200,'r','.')
xlabel('fast LM activity in fb animals - slope');
ylabel('fb-com-slope');

p_ov=signrank([allanimalsLMgo(1:6),allanimalsLMnogo(1:6)],...
    [getallvals(fbreduced.animaltau_com_go,op,type,deltat),getallvals(fbreduced.animaltau_com_nogo,op,type,deltat)]);
hold on; text(-6,0,['pval signrank = ',num2str(p_ov)])

%% pre 200 and post 200 silencing onsets first 3 (after 0) or last 3 : first 3: 65 120 190, last 3:260,320,390

op = 'slope'; % slope, norml1, norml8
deltat = 0.065;
propname = op;

if strcmp(op,'slope')
yl = [-10,0];
else
yl = [-0.5,1.5];
end


%  same as above but in 2 columns
figure;
type = 'first_perlag';
hold on;multiplebars(getallvals(fb150.animaltau_com_go,op,type,deltat),...
    getallvals(fb150.animaltau_com_nogo,op,type,deltat),1,2,yl,'g','r')

type = 'last_perlag';
hold on;multiplebars(getallvals(fb150.animaltau_com_go,op,type,deltat),...
    getallvals(fb150.animaltau_com_nogo,op,type,deltat),3,4,yl,'g','r')


xticks([1 2 3 4])
xticklabels({'fbcom-go-first','fbcom-nogo-first','fbcom-go-last','fbcom-nogo-last'})

ylabel(propname)
xlim([0 6])
ylim([yl(1),yl(2)+0.6])

%%
function plotgonogopairs(gotarget,gocom,nogotarget,nogocom,...
    ci_gotarget,ci_gocom,ci_nogotarget,ci_nogocom)

scatter(gotarget,gocom,100,'g.');hold on;scatter(nogotarget,nogocom,100,'r.');hold on;

errorbar(gotarget,gocom,ci_gotarget,'horizontal','g','LineStyle','none');hold on;
errorbar(nogotarget,nogocom,ci_nogotarget,'horizontal','r','LineStyle','none');hold on;

errorbar(gotarget,gocom,ci_gocom,'g','LineStyle','none');hold on;
errorbar(nogotarget,nogocom,ci_nogocom,'r','LineStyle','none');hold on;

ax=gca;
set(ax, 'xscale','log')
set(ax, 'yscale','log')
xl = get(ax,'XLim');
yl = get(ax,'YLim');
hold on;line(ax,[0.1,max(xl(2),yl(2))],[0.1,max(xl(2),yl(2))])
xlabel('acticity-target'); ylabel('com');
end

%%
function multiplelines(a,b,x1,x2,yl)
for i = 1:length(a)
    hold on;line([x1,x2],[a(i),b(i)],'color','k')
    hold on; scatter([x1,x2],[a(i),b(i)],100,'k.')
end

pff = signrank(a,b);
hold on;text(x1,yl(2)+.05*range(yl),['signrank p = ',num2str(pff)])
pff = ranksum(a,b);
hold on;text(x1,yl(2),['ranksum p = ',num2str(pff)])
end
%%
function u_multiplebars(u,g,n,xu,yl,c)
hold on;
hold on;bar(xu,nanmean(u),0.5,c,'EdgeColor','none')
hold on;errorbar(xu,nanmean(u),2*nanstd(u)/sqrt(numel(find(~isnan(u)))),'k')

pg = ranksum(u,g);
pn = ranksum(u,n);
hold on;text(xu,yl(2),['ranksum p go = ',num2str(pg)])
hold on;text(xu,yl(2)+.05*range(yl),['ranksum p nogo = ',num2str(pn)])

end
%%
function multiplebars(a,b,x1,x2,yl,c1,c2)
hold on;bar(x1,nanmean(a),0.5,c1,'EdgeColor','none')
hold on;errorbar(x1,nanmean(a),2*nanstd(a)/sqrt(numel(find(~isnan(a)))),'k')
hold on;bar(x2,nanmean(b),0.5,c2,'EdgeColor','none')
hold on;errorbar(x2,nanmean(b),2*nanstd(b)/sqrt(numel(find(~isnan(b)))),'k')

pff = signrank(a,b);
hold on;text(x1,yl(2)+.05*range(yl),['signrank p = ',num2str(pff)])
pff = ranksum(a,b);
hold on;text(x1,yl(2),['ranksum p = ',num2str(pff)])
end
%%
function out = getallvals(st,op,type,deltat) % op: 'slope','norml1','norml8',type: 'peranimal','perlag'
% for peranimal:
% lag0 shoud have 8 elements, and no nans, lag1 should have 7 non-nan
% values and lag8, 1. So, lag1 and lag 8 are nanmeaned, but lag0 is meaned.
% meaning that if less than 8 diagonal value in the animal, that animal
% would return nan and would not be included
if strcmp(op,'norml1')
    if strcmp(type,'peranimal')
        out = cellfun(@(x) nanmean(x.lag1)/mean(x.lag0), st);
    elseif strcmp(type,'perlag')
        out = cell2mat(cellfun(@(x) x.lag1/x.lag0, st,'UniformOutput',0)')';
    end
elseif strcmp(op,'norml8')
    if strcmp(type,'peranimal')
        out = cellfun(@(x) nanmean(x.lag8)/mean(x.lag0), st);
    elseif strcmp(type,'perlag')
        out = cell2mat(cellfun(@(x) x.lag8/x.lag0, st,'UniformOutput',0)')';
    end
elseif strcmp(op,'slope')
    if strcmp(type,'peranimal')
        out = cellfun(@(x) (nanmean(x.lag1)-mean(x.lag0))/deltat, st);
    elseif strcmp(type,'perlag')
        out = cell2mat(cellfun(@(x) (x.lag1-x.lag0)/deltat, st,'UniformOutput',0)')';
    elseif strcmp(type,'first_perlag')
        out = cell2mat(cellfun(@(x) (x.lag1(2:4)-x.lag0(2:4))/deltat, st,'UniformOutput',0)')';
    elseif strcmp(type,'last_perlag')
        out = cell2mat(cellfun(@(x) (x.lag1(5:7)-x.lag0(5:7))/deltat, st,'UniformOutput',0)')';       
    end
end
end

