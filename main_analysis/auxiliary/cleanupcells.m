% this function takes the cells from the loaded mat file, right after
% loading, and returns them cleaned up based on the specified critera
% In this vesrion, the same criteria is used for V1 and LM, but the input
% can be given as an array
% cleaning up entails putting some fields to nan, not removing any cells.
% If a certain laser bin is cleaned up due to low activity for
% example, all its analysis fields turn to nan.

function [LMcells,V1cells,Behcells,params] = cleanupcells(cleanupcriteria,LMcells,V1cells,Behcells,params)

targetcellnamelist = {'V1cells','LMcells'};

if cleanupcriteria.doclean
    
    %%%%%%%%%%%%%%%%% min activity in a bin to quantify silencing - smb
    if ~isnan(cleanupcriteria.smb_activity)
        % for V1 and LM
        for areaind = 1:2
            
            targetcellname = targetcellnamelist{areaind};
            targetcell = eval(targetcellname);
            % threshold on baseline activity in the bin
            TH = cleanupcriteria.smb_activity*params.edgestep_an/30/1000;
            for i=1:length(targetcell)
                % go
                for l=1:length(targetcell{i}.laAbs.go)
                    if nanmean(targetcell{i}.laAbs.go{l}) <= TH
                        targetcell{i}.smb.go(l) = nan;
                        if isfield(targetcell{i},'smb_ctrl')
                            targetcell{i}.smb_ctrl.go(l) = nan;
                        end
                        targetcell{i}.smb_p.go(l) = nan;
                        targetcell{i}.smb_z.go(l) = nan;
                        targetcell{i}.smb_bs.go(l) = nan;
                        targetcell{i}.smb_ls.go(l) = nan;
                        targetcell{i}.smb_delta.go(l) = nan;
                        
                        targetcell{i}.laAbs.go{l} = nan(size(targetcell{i}.laAbs.go{l}));
                        targetcell{i}.laAls.go{l} = nan(size(targetcell{i}.laAls.go{l}));
                        if isfield(targetcell{i},'laAbs_ctrl')  && isfield(targetcell{i},'laAls_ctrl')
                            targetcell{i}.laAbs_ctrl.go{l} = nan(size(targetcell{i}.laAbs_ctrl.go{l}));
                            targetcell{i}.laAls_ctrl.go{l} = nan(size(targetcell{i}.laAls_ctrl.go{l}));
                        end
                    end
                end
                % nogo
                for l=1:length(targetcell{i}.laAbs.nogo)
                    if nanmean(targetcell{i}.laAbs.nogo{l}) <= TH
                        targetcell{i}.smb.nogo(l) = nan;
                       if isfield(targetcell{i},'smb_ctrl')
                            targetcell{i}.smb_ctrl.nogo(l) = nan;
                        end
                        targetcell{i}.smb_p.nogo(l) = nan;
                        targetcell{i}.smb_z.nogo(l) = nan;
                        targetcell{i}.smb_bs.nogo(l) = nan;
                        targetcell{i}.smb_ls.nogo(l) = nan;
                        targetcell{i}.smb_delta.nogo(l) = nan;
                        
                        targetcell{i}.laAbs.nogo{l} = nan(size(targetcell{i}.laAbs.nogo{l}));
                        targetcell{i}.laAls.nogo{l} = nan(size(targetcell{i}.laAls.nogo{l}));
                        if isfield(targetcell{i},'laAbs_ctrl')  && isfield(targetcell{i},'laAls_ctrl')
                            targetcell{i}.laAbs_ctrl.nogo{l} = nan(size(targetcell{i}.laAbs_ctrl.nogo{l}));
                            targetcell{i}.laAls_ctrl.nogo{l} = nan(size(targetcell{i}.laAls_ctrl.nogo{l}));
                        end
                    end
                end
            end
            % put targetcell back into V1cells or LMcells
            eval(sprintf('%s = targetcell;',targetcellname));
        end
    end
    %%%%%%%%%%%%%%%%%
    if ~isnan(cleanupcriteria.significance)
    end
    
    
end