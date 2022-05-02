function angs = chancecor(mat)
% n neurons * 8
% mat = cell2mat(animalmodel{1}.lmodel.rep{1}.part{1}.go.vec_n);
%mat = cell2mat(animalmodel{1}.lmodel.rep{1}.part{1}.go.vec_n);
allvals = reshape(mat,1,[]);

nrep=1000;
randvec = nan(size(mat,1),nrep);
for n = 1:size(mat,1)
   for rep = 1:nrep
       %randvec(n,rep) = allvals(randsample(length(allvals),1));
       %randvec(n,rep) = (rand(1)*(max(mat(n,:))-min(mat(n,:)))) + min(mat(n,:));
       randvec(n,rep) = rand(1);
   end
end


nrep =1000;
angs = nan(1,nrep);
for rep = 1:nrep
    ind = randsample(size(randvec,2),2);
    angs(rep) =  (randvec(:,ind(1))/norm(randvec(:,ind(1))))' * ...
        (randvec(:,ind(2))/norm(randvec(:,ind(2))));

end
%figure;hist(angs)


