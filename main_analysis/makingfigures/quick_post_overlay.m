
% root = '';

fb = ...
    openfig([root,'FF150_tc_go.fig']);

fbs = ...
    openfig([root,'simul_ff_tc_go.fig']);
%%%
% blue fast

l=findobj(fbs,'type','Line');
l.Color = [0.5 0.5 0.5];
l=findobj(fbs,'type','Patch');
l.FaceColor = [0.5 0.5 0.5];
l=findobj(fbs,'type','Errorbar');
l.Color = [0.5 0.5 0.5];
l=findobj(fbs,'type','Text');
l.Color = [0.5 0.5 0.5];

% combine

copyobj(findobj(fbs,'type','Line'),findobj(fb,'type','axes'))
copyobj(findobj(fbs,'type','Patch'),findobj(fb,'type','axes'))
copyobj(findobj(fbs,'type','Errorbar'),findobj(fb,'type','axes'))
copyobj(findobj(fbs,'type','Text'),findobj(fb,'type','axes'))
