function inds = findnonnancells_go(animalmodel,animalnum,rep,p)
% mixed go and nogo
out = [];
for type = 0:1
lag = 1;
a1=sum(animalmodel{animalnum}.lmodel.rep{rep}.part{p}.go.X{lag}...
    (find(animalmodel{animalnum}.lmodel.rep{rep}.part{p}.go.Y{lag}==type),:),1);
lag=2;
a2=sum(animalmodel{animalnum}.lmodel.rep{rep}.part{p}.go.X{lag}...
    (find(animalmodel{animalnum}.lmodel.rep{rep}.part{p}.go.Y{lag}==type),:),1);
lag=3;
a3=sum(animalmodel{animalnum}.lmodel.rep{rep}.part{p}.go.X{lag}...
    (find(animalmodel{animalnum}.lmodel.rep{rep}.part{p}.go.Y{lag}==type),:),1);
lag=4;
a4=sum(animalmodel{animalnum}.lmodel.rep{rep}.part{p}.go.X{lag}...
    (find(animalmodel{animalnum}.lmodel.rep{rep}.part{p}.go.Y{lag}==type),:),1);
lag=5;
a5=sum(animalmodel{animalnum}.lmodel.rep{rep}.part{p}.go.X{lag}...
    (find(animalmodel{animalnum}.lmodel.rep{rep}.part{p}.go.Y{lag}==type),:),1);
lag=6;
a6=sum(animalmodel{animalnum}.lmodel.rep{rep}.part{p}.go.X{lag}...
    (find(animalmodel{animalnum}.lmodel.rep{rep}.part{p}.go.Y{lag}==type),:),1);
lag=7;
a7=sum(animalmodel{animalnum}.lmodel.rep{rep}.part{p}.go.X{lag}...
    (find(animalmodel{animalnum}.lmodel.rep{rep}.part{p}.go.Y{lag}==type),:),1);
lag=8;
a8=sum(animalmodel{animalnum}.lmodel.rep{rep}.part{p}.go.X{lag}...
    (find(animalmodel{animalnum}.lmodel.rep{rep}.part{p}.go.Y{lag}==type),:),1);

out=[out;sum([a1;a2;a3;a4;a5;a6;a7;a8],1)];

end

inds = find(~isnan(sum(out,1)));