% this function takes the cells from the loaded mat file, right after
% loading, and removes some cells based on the specified criteria, the same
% criteria is used for V1 and LM, but the input can be given as an array
% if cell_keep_ind is given as a non empty struct, the function just
% indexes LM and V1 cells based on cell_keep_ind, and does not recalculate
% any indexes.

% cell_keep_ind is a struct with V1cells and LMcells fields

function [cell_keep_ind,LMcells,V1cells,Behcells,params] = removecells(cell_keep_ind,cellremovecriteria,LMcells,V1cells,Behcells,params)
%% parans
n_example_cells = 10;
%%
targetcellnamelist = {'V1cells','LMcells'};

if cellremovecriteria.doclean
    
    % idindexes are given: just index the cells
    if  isstruct(cell_keep_ind)%~(isempty(cell_keep_ind) || isnan( cell_keep_ind))
        
        for areaind = 1:2
            
            targetcellname = targetcellnamelist{areaind};
            targetcell = eval(targetcellname);
            
            % index it
            targetcell = targetcell(getfield(cell_keep_ind,targetcellname));
            
            % put targetcell back into V1cells or LMcells
            eval(sprintf('%s = targetcell;',targetcellname));
            
        end
        
        
        % if index are not given (empty or nan) - get indexes
    else
        % final indices, filtered through all criteria: The order of
        % criteria doesn't matter, outcome is the and of all.
        cell_keep_ind = struct;
        cell_keep_ind.V1cells = 1:1:length(V1cells);
        cell_keep_ind.LMcells = 1:1:length(LMcells);
        % if removing based on responsiveness
        if ~isnan(cellremovecriteria.responsiveness)
            % for each area, put it in targetcells, filter
            for areaind = 1:2
                clear keep_ind
                targetcellname = targetcellnamelist{areaind};
                targetcell = eval(targetcellname);
                %%% get resp measure per go and nogo separately
                if cellremovecriteria.responsiveness(1) == 1 % if method is 1
                    nogoind=find(cellfun(@(x) signrank(x.stimspikes.nogo,x.prestimspikes.nogo)<=cellremovecriteria.sigTH , targetcell));
                    goind=find(cellfun(@(x) signrank(x.stimspikes.go,x.prestimspikes.go)<=cellremovecriteria.sigTH , targetcell));
                elseif cellremovecriteria.responsiveness(1) == 2 % if method is 2
                    nogoind=find(cellfun(@(x) nanmean(x.stimspikes.nogo-x.prestimspikes.nogo)>0 , targetcell));
                    goind=find(cellfun(@(x) nanmean(x.stimspikes.go-x.prestimspikes.go)>0 , targetcell));                   
                elseif cellremovecriteria.responsiveness(1) == 3 % if method is 3
                    nogoind=intersect(find(cellfun(@(x) nanmean(x.stimspikes.nogo-x.prestimspikes.nogo)>0 , targetcell)),...
                        find(cellfun(@(x) signrank(x.stimspikes.nogo,x.prestimspikes.nogo)<=cellremovecriteria.sigTH , targetcell)));
                    goind=intersect(find(cellfun(@(x) nanmean(x.stimspikes.go-x.prestimspikes.go)>0 , targetcell)),...
                        find(cellfun(@(x) signrank(x.stimspikes.go,x.prestimspikes.go)<=cellremovecriteria.sigTH , targetcell)));                   
               elseif cellremovecriteria.responsiveness(1) == 4 % if method is 3
                    nogoind=intersect(find(cellfun(@(x) nanmean(x.stimspikes.nogo)>cellremovecriteria.minspkTH , targetcell)),...
                        find(cellfun(@(x) signrank(x.stimspikes.nogo,x.prestimspikes.nogo)<=cellremovecriteria.sigTH , targetcell)));
                    goind=intersect(find(cellfun(@(x) nanmean(x.stimspikes.go)>cellremovecriteria.minspkTH , targetcell)),...
                        find(cellfun(@(x) signrank(x.stimspikes.go,x.prestimspikes.go)<=cellremovecriteria.sigTH , targetcell)));                   
                else
                    error('Unknown method')
                end
                %%% combine go and nogo rep mesures
                if cellremovecriteria.responsiveness(2) == 1  
                    keep_ind = union(nogoind,goind);
                elseif cellremovecriteria.responsiveness(2) == 2 
                    keep_ind = intersect(nogoind,goind);
                else
                    error('Unknown method')
                end
                
                % intersect the cell_keep_ind field with the new keep_ind
                % from this section
                cell_keep_ind = setfield(cell_keep_ind,targetcellname,intersect(keep_ind,getfield(cell_keep_ind,targetcellname)));                
            end
        end
            
        if ~isnan(cellremovecriteria.stability)
            % for each area, put it in targetcells, filter
            for areaind = 1:2
                clear keep_ind
                targetcellname = targetcellnamelist{areaind};
                targetcell = eval(targetcellname);  
                
                lindrift_go=(cellfun(@(x) corrcoef(x.stimspikes.go,1:length(x.stimspikes.go)) , targetcell,'UniformOutput',0));
                lindrift_nogo=(cellfun(@(x) corrcoef(x.stimspikes.nogo,1:length(x.stimspikes.nogo)) , targetcell,'UniformOutput',0));
                % average of linear drift in go and nogo trials is used as
                % a measure of stability:
                lindrift_all = cellfun(@(x,y) (x(1,2)+y(1,2))/2, lindrift_go,lindrift_nogo);
                %
                keep_ind = find( abs(lindrift_all) < cellremovecriteria.lindriftTH );
                
                % intersect the cell_keep_ind field with the new keep_ind
                % from this section
                cell_keep_ind = setfield(cell_keep_ind,targetcellname,intersect(keep_ind,getfield(cell_keep_ind,targetcellname)));                
            end
                             
        end
        
        
        %%%%%%%%%%%% after finishing all criteria:
        %%%%%%%%%%%% modify targetcell and put targetcell back into V1cells or LMcells
        % plot stuff
        for areaind = 1:2
            targetcellname = targetcellnamelist{areaind};
            targetcell = eval(targetcellname);
            

            not_targetcell = targetcell(setdiff(1:1:length(targetcell),getfield(cell_keep_ind,targetcellname)));            
            targetcell = targetcell(getfield(cell_keep_ind,targetcellname));
            
            if cellremovecriteria.showexampleplots
                figure;
                for ex_c = 1:n_example_cells
                    targetex = randi(size(targetcell,2));
                    subplot(n_example_cells,4,4*(ex_c-1)+1);                
                    imagesc(targetcell{targetex}.nbs.go);
                    subplot(n_example_cells,4,4*(ex_c-1)+2);
                    imagesc(targetcell{targetex}.nbs.nogo);
                    
                    not_targetex = randi(size(not_targetcell,2));
                    subplot(n_example_cells,4,4*(ex_c-1)+3);
                    imagesc(not_targetcell{not_targetex}.nbs.go);
                    subplot(n_example_cells,4,4*(ex_c-1)+4);
                    imagesc(not_targetcell{not_targetex}.nbs.nogo);
                end
                
                
                
            end
            
            eval(sprintf('%s = targetcell;',targetcellname));
            fprintf('%d cells were removed out of %d in %s\n',...
                size(not_targetcell,2),size(not_targetcell,2)+size(targetcell,2),targetcellnamelist{areaind} )
        end
    end
end