function [Cmatgo_mu,Cmatnogo_mu,Crepgo_mu,Crepnogo_mu,C_go_mu,C_nogo_mu] = X_prepareCmatrices_mu(animalmodel,includedanimals)
nlags = 8;

Cmatgo_mu = nan(8,8,includedanimals(end));
Cmatnogo_mu = nan(8,8,includedanimals(end));
Crepgo_mu = nan(8,8,animalmodel{1}.lmodel.nrep);
Crepnogo_mu = nan(8,8,animalmodel{1}.lmodel.nrep);

% for each animal, average of reps is calculated
for animalnum = includedanimals
    for repi = 1: animalmodel{animalnum}.lmodel.nrep
        for i = 1:8
            for j = 1:8
                % comvec part1 (i) * comvec part2 (j)
            
                % mu vecs
                Crepgo_mu(i,j,repi) = ...
                    (animalmodel{animalnum}.lmodel.rep{repi}.part{1}.go.mu_dist{i}' * animalmodel{animalnum}.lmodel.rep{repi}.part{2}.go.mu_dist{j} + ...
                    animalmodel{animalnum}.lmodel.rep{repi}.part{2}.go.mu_dist{i}' * animalmodel{animalnum}.lmodel.rep{repi}.part{1}.go.mu_dist{j})/2;
                Crepnogo_mu(i,j,repi) = ...
                    (animalmodel{animalnum}.lmodel.rep{repi}.part{1}.nogo.mu_dist{i}' * animalmodel{animalnum}.lmodel.rep{repi}.part{2}.nogo.mu_dist{j} + ...
                    animalmodel{animalnum}.lmodel.rep{repi}.part{2}.nogo.mu_dist{i}' * animalmodel{animalnum}.lmodel.rep{repi}.part{1}.nogo.mu_dist{j})/2;
           
            end
        end
    end
    
     Cmatgo_mu(:,:,animalnum) = nanmean(Crepgo_mu,3);
    Cmatnogo_mu(:,:,animalnum) = nanmean(Crepnogo_mu,3);
end

%%%%%%%%%%%%%%%%%%
% (nlags*nanimals)*7 or 8 (lags)
C_go_mu = nan(0,nlags);
C_nogo_mu = nan(0,nlags);
count = 1;
for animalnum = includedanimals
    for i = 1:nlags
        C_go_mu(end+1,:) = nan(1,nlags);
        C_nogo_mu(end+1,:) = nan(1,nlags);
        for j = 1:nlags
            if j>=i
%                 C_go_lda(count,j-i+1) = Cmatgo_lda(i,j,animalnum);
%                 C_nogo_lda(count,j-i+1) = Cmatnogo_lda(i,j,animalnum);
%                 count = count + 1;
                C_go_mu(end,j-i+1) = Cmatgo_mu(i,j,animalnum);
                C_nogo_mu(end,j-i+1) = Cmatnogo_mu(i,j,animalnum);               
            end
        end
    end
end