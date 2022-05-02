
root = '/mnt/data/Mitra/figs/P2_L/bothdircombined/2020_02_24/LDA/version2/cosinesimilaritymats_and_timeconstants/saveres/';

lmslow = ...
    openfig([root,'LM-slow0.8_Activity_ff_fb_combined65Xstyle2nrep100_2021_01_10_22_41_32/LM-slow0.8-raw-65ms-combiend.fig']);

lmfast = ...
    openfig([root,'LM-fast0.2_Activity_ff_fb_combined65Xstyle2nrep100_2021_01_10_22_39_16/LM-fast0.2-raw-65ms-combiend.fig']);

v1slow = ...
    openfig([root,'V1-slow0.8_Activity_ff_fb_combined65Xstyle2nrep100_2021_01_10_22_43_31/V1-slow0.8-raw-65ms-combiend.fig']);

v1fast = ...
    openfig([root,'V1-fast0.2_Activity_ff_fb_combined65Xstyle2nrep100_2021_01_10_22_45_43/V1-fast0.2-raw-65ms-combiend.fig']);
%%
% blue fast

l=findobj(lmfast,'type','Line');
l.Color='c';
l=findobj(lmfast,'type','Patch');
l.FaceColor='c';
l=findobj(lmfast,'type','Errorbar');
l.Color='c';
l=findobj(lmfast,'type','Text');
l.Color='c';


l=findobj(v1fast,'type','Line');
l.Color='c';
l=findobj(v1fast,'type','Patch');
l.FaceColor='c';
l=findobj(v1fast,'type','Errorbar');
l.Color='c';
l=findobj(v1fast,'type','Text');
l.Color='c';

% combine

copyobj(findobj(lmfast,'type','Line'),findobj(lmslow,'type','axes'))
copyobj(findobj(lmfast,'type','Patch'),findobj(lmslow,'type','axes'))
copyobj(findobj(lmfast,'type','Errorbar'),findobj(lmslow,'type','axes'))
copyobj(findobj(lmfast,'type','Text'),findobj(lmslow,'type','axes'))

copyobj(findobj(v1fast,'type','Line'),findobj(v1slow,'type','axes'))
copyobj(findobj(v1fast,'type','Patch'),findobj(v1slow,'type','axes'))
copyobj(findobj(v1fast,'type','Errorbar'),findobj(v1slow,'type','axes'))
copyobj(findobj(v1fast,'type','Text'),findobj(v1slow,'type','axes'))
