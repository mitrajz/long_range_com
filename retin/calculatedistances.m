function [retin_V1,retin_LM] = calculatedistances(retin_V1,retin_LM)
% units: 1 unit in x: size of retinotopy block in X, 1 unit y: size of
% patch in Y
% The way distance is calculated assumes they are the same. To convert
% distance to cm, should multiply by the length of match in cm, and to
% degrees, should do a rough conversion of cm to degrees (1 edge length ~
% 20 deg around the center)
% IMPOTANT: note that the center is center of 2 shanks, so the average
% x or Y each shank would not be zero


V1x = nanmean([cell2mat(retin_V1{1}.retinmap.xgcenter),cell2mat(retin_V1{2}.retinmap.xgcenter)]);
V1y = nanmean([cell2mat(retin_V1{1}.retinmap.ygcenter),cell2mat(retin_V1{2}.retinmap.ygcenter)]);
LMx = nanmean([cell2mat(retin_LM{1}.retinmap.xgcenter),cell2mat(retin_LM{2}.retinmap.xgcenter)]);
LMy = nanmean([cell2mat(retin_LM{1}.retinmap.ygcenter),cell2mat(retin_LM{2}.retinmap.ygcenter)]);

% for V1
for shank = 1:2
    % set local distances
    retin_V1{shank}.retinmap.distancelocal.x = cellfun(@(x) x-V1x,retin_V1{shank}.retinmap.xgcenter,'UniformOutput',0);
    retin_V1{shank}.retinmap.distancelocal.y = cellfun(@(x) x-V1y,retin_V1{shank}.retinmap.ygcenter,'UniformOutput',0);
    retin_V1{shank}.retinmap.distancelocal.distance = cellfun(@(x,y) sqrt(x^2 + y^2), ...
        retin_V1{shank}.retinmap.distancelocal.x,retin_V1{shank}.retinmap.distancelocal.y,'UniformOutput',0);
    % set long range diatances
    retin_V1{shank}.retinmap.distancelong.x = cellfun(@(x) x-LMx,retin_V1{shank}.retinmap.xgcenter,'UniformOutput',0);
    retin_V1{shank}.retinmap.distancelong.y = cellfun(@(x) x-LMy,retin_V1{shank}.retinmap.ygcenter,'UniformOutput',0);
    retin_V1{shank}.retinmap.distancelong.distance = cellfun(@(x,y) sqrt(x^2 + y^2), ...
        retin_V1{shank}.retinmap.distancelong.x,retin_V1{shank}.retinmap.distancelong.y,'UniformOutput',0);
      
end

% for LM
for shank = 1:2
    % set local distances
    retin_LM{shank}.retinmap.distancelocal.x = cellfun(@(x) x-LMx,retin_LM{shank}.retinmap.xgcenter,'UniformOutput',0);
    retin_LM{shank}.retinmap.distancelocal.y = cellfun(@(x) x-LMy,retin_LM{shank}.retinmap.ygcenter,'UniformOutput',0);
    retin_LM{shank}.retinmap.distancelocal.distance = cellfun(@(x,y) sqrt(x^2 + y^2), ...
        retin_LM{shank}.retinmap.distancelocal.x,retin_LM{shank}.retinmap.distancelocal.y,'UniformOutput',0);
    % set long range diatances
    retin_LM{shank}.retinmap.distancelong.x = cellfun(@(x) x-V1x,retin_LM{shank}.retinmap.xgcenter,'UniformOutput',0);
    retin_LM{shank}.retinmap.distancelong.y = cellfun(@(x) x-V1y,retin_LM{shank}.retinmap.ygcenter,'UniformOutput',0);
    retin_LM{shank}.retinmap.distancelong.distance = cellfun(@(x,y) sqrt(x^2 + y^2), ...
        retin_LM{shank}.retinmap.distancelong.x,retin_LM{shank}.retinmap.distancelong.y,'UniformOutput',0);
      
end