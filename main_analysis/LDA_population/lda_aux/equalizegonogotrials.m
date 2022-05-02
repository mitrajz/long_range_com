function targetcell = equalizegonogotrials(targetcell,frnormalize,cutflag,randeq);
% randeq = 0 : the first n trials are used when equalizing
% randeq = 1: random selection of trials
samelsbstrials = 0;

rng(1,'twister')

lsntrials = nan(1,8);
bsntrials = nan(1,8);
for lag=1:8
    if cutflag
        param=5/4;
    else
        param=1;
    end
    lsntrials(lag) = round(min(min(unique(cellfun(@(x) size(x.laAls.go{lag},1),targetcell))),...
        min(unique(cellfun(@(x) size(x.laAls.nogo{lag},1),targetcell))))/...
        param);
   
    bsntrials(lag) = round(min(min(unique(cellfun(@(x) size(x.laAbs.go{lag},1),targetcell))),...
        min(unique(cellfun(@(x) size(x.laAbs.nogo{lag},1),targetcell))))/...
        param);
    if samelsbstrials
        minntrials = min(lsntrials(lag),bsntrials(lag));
        lsntrials(lag) = minntrials;
        bsntrials(lag) = minntrials;
    end
end



for n= 1:length(targetcell)

    gobs = targetcell{n}.nbs.go;
    nogobs = targetcell{n}.nbs.nogo;
    targetcell{n}.nbs.go = cell(1,8);
    targetcell{n}.nbs.nogo = cell(1,8);
    
    for lag=1:8
        targetcell{n}.nbs.go{lag} = gobs;
        targetcell{n}.nbs.nogo{lag} = nogobs;
        
        if randeq
            % to make sure the same for all cells of the same animals
            rng(1,'twister');
            lasergoind = randperm(length(targetcell{n}.laAls.go{lag}),lsntrials(lag));
            rng(1,'twister');
            lasernogoind = randperm(length(targetcell{n}.laAls.nogo{lag}),lsntrials(lag));
            rng(1,'twister');
            baselinegoind = randperm(length(targetcell{n}.laAbs.go{lag}),bsntrials(lag));
            rng(1,'twister');
            baselinenogoind = randperm(length(targetcell{n}.laAbs.nogo{lag}),bsntrials(lag));
        else
            lasergoind = (length(targetcell{n}.laAls.go{lag})-lsntrials(lag)+1):length(targetcell{n}.laAls.go{lag});
            lasernogoind = (length(targetcell{n}.laAls.nogo{lag})-lsntrials(lag)+1):length(targetcell{n}.laAls.nogo{lag});
            baselinegoind = (length(targetcell{n}.laAbs.go{lag})-bsntrials(lag)+1):length(targetcell{n}.laAbs.go{lag});
            baselinenogoind = (length(targetcell{n}.laAbs.nogo{lag})-bsntrials(lag)+1):length(targetcell{n}.laAbs.nogo{lag});

        end
        % laser go
        targetcell{n}.laAls.go{lag} = targetcell{n}.laAls.go{lag}(lasergoind);
        targetcell{n}.nls.go{lag} = targetcell{n}.nls.go{lag}(lasergoind,:);
        % laser nogo
        targetcell{n}.laAls.nogo{lag} = targetcell{n}.laAls.nogo{lag}(lasernogoind);
        targetcell{n}.nls.nogo{lag} = targetcell{n}.nls.nogo{lag}(lasernogoind,:);
        % baseline go
        targetcell{n}.laAbs.go{lag} = targetcell{n}.laAbs.go{lag}(baselinegoind);
        targetcell{n}.nbs.go{lag} = targetcell{n}.nbs.go{lag}(baselinegoind,:);
        % baseline nogo
        targetcell{n}.laAbs.nogo{lag} = targetcell{n}.laAbs.nogo{lag}(baselinenogoind);
        targetcell{n}.nbs.nogo{lag} = targetcell{n}.nbs.nogo{lag}(baselinenogoind,:);
        
          
        if frnormalize
            targetcell{n}.laAls.go{lag} = (targetcell{n}.laAls.go{lag}-targetcell{n}.Avfr);%/targetcell{n}.Avfr;
            targetcell{n}.laAls.nogo{lag} = (targetcell{n}.laAls.nogo{lag}-targetcell{n}.Avfr);%/targetcell{n}.Avfr;
            targetcell{n}.laAbs.go{lag} = (targetcell{n}.laAbs.go{lag}-targetcell{n}.Avfr);%/targetcell{n}.Avfr;
            targetcell{n}.laAbs.nogo{lag} = (targetcell{n}.laAbs.nogo{lag}-targetcell{n}.Avfr);%/targetcell{n}.Avfr;
        end
        
        
    end
end