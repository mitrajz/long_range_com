function targetcell = virtualgonocomb(targetcell)
old_L = length(targetcell);
targetcell = [targetcell targetcell];
for n = 1:old_L
    clear g ng
    g = targetcell{old_L+n}.stimspikes.go;
    ng = targetcell{old_L+n}.stimspikes.nogo;
    targetcell{old_L+n}.stimspikes.go = ng;
    targetcell{old_L+n}.stimspikes.nogo = g;
    
     clear g ng
    g = targetcell{old_L+n}.prestimspikes.go;
    ng = targetcell{old_L+n}.prestimspikes.nogo;
    targetcell{old_L+n}.prestimspikes.go = ng;
    targetcell{old_L+n}.prestimspikes.nogo = g;
    
    clear g ng
    g = targetcell{old_L+n}.nbs.go;
    ng = targetcell{old_L+n}.nbs.nogo;
    targetcell{old_L+n}.nbs.go = ng;
    targetcell{old_L+n}.nbs.nogo = g;
    
    clear g ng
    g = targetcell{old_L+n}.nls.go;
    ng = targetcell{old_L+n}.nls.nogo;
    targetcell{old_L+n}.nls.go = ng;
    targetcell{old_L+n}.nls.nogo = g;
    
    clear g ng
    g = targetcell{old_L+n}.nbsAv.go;
    ng = targetcell{old_L+n}.nbsAv.nogo;
    targetcell{old_L+n}.nbsAv.go = ng;
    targetcell{old_L+n}.nbsAv.nogo = g;
    
    clear g ng
    g = targetcell{old_L+n}.nlsAv.go;
    ng = targetcell{old_L+n}.nlsAv.nogo;
    targetcell{old_L+n}.nlsAv.go = ng;
    targetcell{old_L+n}.nlsAv.nogo = g;
    
    clear g ng
    g = targetcell{old_L+n}.smb.go;
    ng = targetcell{old_L+n}.smb.nogo;
    targetcell{old_L+n}.smb.go = ng;
    targetcell{old_L+n}.smb.nogo = g;
    
    clear g ng
    g = targetcell{old_L+n}.smb_p.go;
    ng = targetcell{old_L+n}.smb_p.nogo;
    targetcell{old_L+n}.smb_p.go = ng;
    targetcell{old_L+n}.smb_p.nogo = g;
    
    clear g ng
    g = targetcell{old_L+n}.smb_z.go;
    ng = targetcell{old_L+n}.smb_z.nogo;
    targetcell{old_L+n}.smb_z.go = ng;
    targetcell{old_L+n}.smb_z.nogo = g;
    
    clear g ng
    g = targetcell{old_L+n}.smb_bs.go;
    ng = targetcell{old_L+n}.smb_bs.nogo;
    targetcell{old_L+n}.smb_bs.go = ng;
    targetcell{old_L+n}.smb_bs.nogo = g;
    
    clear g ng
    g = targetcell{old_L+n}.smb_ls.go;
    ng = targetcell{old_L+n}.smb_ls.nogo;
    targetcell{old_L+n}.smb_ls.go = ng;
    targetcell{old_L+n}.smb_ls.nogo = g;
    
    clear g ng
    g = targetcell{old_L+n}.smb_delta.go;
    ng = targetcell{old_L+n}.smb_delta.nogo;
    targetcell{old_L+n}.smb_delta.go = ng;
    targetcell{old_L+n}.smb_delta.nogo = g;
    
    clear g ng
    g = targetcell{old_L+n}.smb_centers.go;
    ng = targetcell{old_L+n}.smb_centers.nogo;
    targetcell{old_L+n}.smb_centers.go = ng;
    targetcell{old_L+n}.smb_centers.nogo = g;
    
    clear g ng
    g = targetcell{old_L+n}.laAbs.go;
    ng = targetcell{old_L+n}.laAbs.nogo;
    targetcell{old_L+n}.laAbs.go = ng;
    targetcell{old_L+n}.laAbs.nogo = g;
    
    clear g ng
    g = targetcell{old_L+n}.laAls.go;
    ng = targetcell{old_L+n}.laAls.nogo;
    targetcell{old_L+n}.laAls.go = ng;
    targetcell{old_L+n}.laAls.nogo = g;
    
    clear g ng
    g = targetcell{old_L+n}.effectsize_lags.go;
    ng = targetcell{old_L+n}.effectsize_lags.nogo;
    targetcell{old_L+n}.effectsize_lags.go = ng;
    targetcell{old_L+n}.effectsize_lags.nogo = g;
    
    clear g ng
    g = targetcell{old_L+n}.effectsize_lags_pre.go;
    ng = targetcell{old_L+n}.effectsize_lags_pre.nogo;
    targetcell{old_L+n}.effectsize_lags_pre.go = ng;
    targetcell{old_L+n}.effectsize_lags_pre.nogo = g;
    
    clear g ng
    g = targetcell{old_L+n}.pval_lags.go;
    ng = targetcell{old_L+n}.pval_lags.nogo;
    targetcell{old_L+n}.pval_lags.go = ng;
    targetcell{old_L+n}.pval_lags.nogo = g;
    
end