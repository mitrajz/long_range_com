function [targetcell] = subtract_prestim(targetcell,binsizems)

for i=1:length(targetcell)
    pre = str2num(binsizems)*(nanmean([targetcell{i}.prestimspikes.go;targetcell{i}.prestimspikes.nogo])/500);
    % go
    for l=1:length(targetcell{i}.laAbs.go)
        targetcell{i}.laAbs.go{l} = targetcell{i}.laAbs.go{l} - pre;
        targetcell{i}.laAls.go{l} = targetcell{i}.laAls.go{l} - pre;
        
        targetcell{i}.smb.go(l) = targetcell{i}.smb.go(l) - pre;
        targetcell{i}.smb_bs.go(l) = targetcell{i}.smb_bs.go(l) - pre;
        targetcell{i}.smb_ls.go(l) = targetcell{i}.smb_ls.go(l) - pre;
       
    end
    % nogo
    for l=1:length(targetcell{i}.laAbs.go)
        targetcell{i}.laAbs.nogo{l} = targetcell{i}.laAbs.nogo{l} - pre;
        targetcell{i}.laAls.nogo{l} = targetcell{i}.laAls.nogo{l} - pre;
        
        targetcell{i}.smb.nogo(l) = targetcell{i}.smb.nogo(l) - pre;
        targetcell{i}.smb_bs.nogo(l) = targetcell{i}.smb_bs.nogo(l) - pre;
        targetcell{i}.smb_ls.nogo(l) = targetcell{i}.smb_ls.nogo(l) - pre;
    end
end

end