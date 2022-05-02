function [acmatgo,acmatnogo,acrepgo,acrepnogo,ac_go,ac_nogo] = X_prepareCmatrices_notarget_activity(animalmodel,includedanimals,...
    lmodel,V1cells,LMcells,normalize)
% normalize, subtracts the average during-stimulus activity fromeacg
% neurons

nlags = 8;

acmatgo = nan(8,8,includedanimals(end));
acmatnogo = nan(8,8,includedanimals(end));
acrepgo = nan(8,8,animalmodel{1}.lmodel.nrep);
acrepnogo = nan(8,8,animalmodel{1}.lmodel.nrep);

tfactor= str2num(lmodel.binsizems)/1000;


for animalnum = includedanimals
     if strcmp(lmodel.exptype,'FB')
        notargetcell = LMcells(find(cellfun(@(x) x.simulcode,LMcells) == animalnum ));
    elseif strcmp(lmodel.exptype,'FF')
        notargetcell = V1cells(find(cellfun(@(x) x.simulcode,V1cells) == animalnum ));
    end
    
    if lmodel.eqtrials
        randeq = 0; % if zero, takes last. Seed is fixed  for random sampling
        notargetcell = equalizegonogotrials(notargetcell,0,0,randeq);% frnormalize
    end
    
    % this 2 vectors have the length of the number of cells. They are the
    % average firing rates of cells during stimulus for all cells (normalize to 0.5s of stim to convert to Hz)
    goAv = cellfun(@(x) nanmean(x.stimspikes.go),notargetcell)/0.5;
    nogoAv = cellfun(@(x) nanmean(x.stimspikes.nogo),notargetcell)/0.5;
    
    for repi = 1: animalmodel{animalnum}.lmodel.nrep
        for i = 1:8
            for j = 1:8
                % activity vec part1 (i) * activiy vec part2 (j)
                
                % each of the 4 vectors are n dim (num of cells), and can
                % contain nans.
                % all activityvecs are converted to Hz.
                
                %%%%% go: 
                
                vec1 = nanmean(animalmodel{animalnum}.lmodel.rep{repi}.part{1}.go.Xn{i}(find(animalmodel{animalnum}.lmodel.rep{repi}.part{1}.go.Yn{i}==0),:),1)/tfactor;
                vec2 = nanmean(animalmodel{animalnum}.lmodel.rep{repi}.part{2}.go.Xn{j}(find(animalmodel{animalnum}.lmodel.rep{repi}.part{2}.go.Yn{j}==0),:),1)/tfactor;
                vec3 = nanmean(animalmodel{animalnum}.lmodel.rep{repi}.part{2}.go.Xn{i}(find(animalmodel{animalnum}.lmodel.rep{repi}.part{2}.go.Yn{i}==0),:),1)/tfactor;
                vec4 = nanmean(animalmodel{animalnum}.lmodel.rep{repi}.part{1}.go.Xn{j}(find(animalmodel{animalnum}.lmodel.rep{repi}.part{1}.go.Yn{j}==0),:),1)/tfactor;
                
                if normalize
                    vec1 = vec1 - goAv;
                    vec2 = vec2 - goAv;
                    vec3 = vec3 - goAv;
                    vec4 = vec4 - goAv;
                    
                end
                
                vec1(find(isnan(vec1))) = 0;
                vec2(find(isnan(vec2))) = 0;
                vec3(find(isnan(vec3))) = 0;
                vec4(find(isnan(vec4))) = 0;
                % set length to 1 
                vec1 = vec1/norm(vec1);
                vec2 = vec2/norm(vec2);
                vec3 = vec3/norm(vec3);
                vec4 = vec4/norm(vec4);
                               
                acrepgo(i,j,repi) = (vec1 * vec2' + vec3 * vec4')/2;
                
                
                %%%%% nogo
                vec1 = nanmean(animalmodel{animalnum}.lmodel.rep{repi}.part{1}.nogo.Xn{i}(find(animalmodel{animalnum}.lmodel.rep{repi}.part{1}.nogo.Yn{i}==0),:),1)/tfactor;
                vec2 = nanmean(animalmodel{animalnum}.lmodel.rep{repi}.part{2}.nogo.Xn{j}(find(animalmodel{animalnum}.lmodel.rep{repi}.part{2}.nogo.Yn{j}==0),:),1)/tfactor;
                vec3 = nanmean(animalmodel{animalnum}.lmodel.rep{repi}.part{2}.nogo.Xn{i}(find(animalmodel{animalnum}.lmodel.rep{repi}.part{2}.nogo.Yn{i}==0),:),1)/tfactor;
                vec4 = nanmean(animalmodel{animalnum}.lmodel.rep{repi}.part{1}.nogo.Xn{j}(find(animalmodel{animalnum}.lmodel.rep{repi}.part{1}.nogo.Yn{j}==0),:),1)/tfactor;
                
                if normalize
                    vec1 = vec1 - nogoAv;
                    vec2 = vec2 - nogoAv;
                    vec3 = vec3 - nogoAv;
                    vec4 = vec4 - nogoAv;
                    
                end
                

                vec1(find(isnan(vec1))) = 0;
                vec2(find(isnan(vec2))) = 0;
                vec3(find(isnan(vec3))) = 0;
                vec4(find(isnan(vec4))) = 0;
                % set length to 1 
                vec1 = vec1/norm(vec1);
                vec2 = vec2/norm(vec2);
                vec3 = vec3/norm(vec3);
                vec4 = vec4/norm(vec4);
                acrepnogo(i,j,repi) = (vec1 * vec2' + vec3 * vec4')/2;
            end
        end
    end
    
    acmatgo(:,:,animalnum) = nanmean(acrepgo,3);
    acmatnogo(:,:,animalnum) = nanmean(acrepnogo,3);
    
end


%%%%%%%%
% (nlags*nanimals)*7 or 8 (lags)
ac_go = nan(0,nlags);
ac_nogo = nan(0,nlags);

for animalnum = includedanimals
    for i = 1:nlags
        ac_go(end+1,:) = nan(1,nlags);
        ac_nogo(end+1,:) = nan(1,nlags);
        for j = 1:nlags
            if j>=i
%                 C_go_lda(count,j-i+1) = Cmatgo_lda(i,j,animalnum);
%                 C_nogo_lda(count,j-i+1) = Cmatnogo_lda(i,j,animalnum);
%                 count = count + 1;
                ac_go(end,j-i+1) = acmatgo(i,j,animalnum);
                ac_nogo(end,j-i+1) = acmatnogo(i,j,animalnum);               
            end
        end
    end
end
