function [Cmatgo_lda,Cmatnogo_lda,Crepgo_lda,Crepnogo_lda,C_go_lda,C_nogo_lda] = X_prepareCmatrices_lda(animalmodel,includedanimals)
nlags = 8;

Cmatgo_lda = nan(8,8,includedanimals(end));
Cmatnogo_lda = nan(8,8,includedanimals(end));
Crepgo_lda = nan(8,8,animalmodel{1}.lmodel.nrep);
Crepnogo_lda = nan(8,8,animalmodel{1}.lmodel.nrep);


for animalnum = includedanimals
    for repi = 1: animalmodel{animalnum}.lmodel.nrep
        for i = 1:8
            for j = 1:8
                % comvec part1 (i) * comvec part2 (j)
                % lda vecs
                Crepgo_lda(i,j,repi) = ...
                    (animalmodel{animalnum}.lmodel.rep{repi}.part{1}.go.vec_n{i}' * animalmodel{animalnum}.lmodel.rep{repi}.part{2}.go.vec_n{j} + ...
                    animalmodel{animalnum}.lmodel.rep{repi}.part{2}.go.vec_n{i}' * animalmodel{animalnum}.lmodel.rep{repi}.part{1}.go.vec_n{j})/2;
                Crepnogo_lda(i,j,repi) = ...
                    (animalmodel{animalnum}.lmodel.rep{repi}.part{1}.nogo.vec_n{i}' * animalmodel{animalnum}.lmodel.rep{repi}.part{2}.nogo.vec_n{j} + ...
                    animalmodel{animalnum}.lmodel.rep{repi}.part{2}.nogo.vec_n{i}' * animalmodel{animalnum}.lmodel.rep{repi}.part{1}.nogo.vec_n{j})/2;
            end
        end
    end
    
    Cmatgo_lda(:,:,animalnum) = nanmean(Crepgo_lda,3);
    Cmatnogo_lda(:,:,animalnum) = nanmean(Crepnogo_lda,3);
    
end


%%%%%%%%
% (nlags*nanimals)*7 or 8 (lags)
C_go_lda = nan(0,nlags);
C_nogo_lda = nan(0,nlags);
count = 1;
for animalnum = includedanimals
    for i = 1:nlags
        C_go_lda(end+1,:) = nan(1,nlags);
        C_nogo_lda(end+1,:) = nan(1,nlags);
        for j = 1:nlags
            if j>=i
%                 C_go_lda(count,j-i+1) = Cmatgo_lda(i,j,animalnum);
%                 C_nogo_lda(count,j-i+1) = Cmatnogo_lda(i,j,animalnum);
%                 count = count + 1;
                C_go_lda(end,j-i+1) = Cmatgo_lda(i,j,animalnum);
                C_nogo_lda(end,j-i+1) = Cmatnogo_lda(i,j,animalnum);               
            end
        end
    end
end
