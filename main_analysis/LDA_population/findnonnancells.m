function inds = findnonnancells(animalmodel,animalnum,rep,p)
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


lag = 1;
b1=sum(animalmodel{animalnum}.lmodel.rep{rep}.part{p}.nogo.X{lag}...
    (find(animalmodel{animalnum}.lmodel.rep{rep}.part{p}.nogo.Y{lag}==type),:),1);
lag=2;
b2=sum(animalmodel{animalnum}.lmodel.rep{rep}.part{p}.nogo.X{lag}...
    (find(animalmodel{animalnum}.lmodel.rep{rep}.part{p}.nogo.Y{lag}==type),:),1);
lag=3;
b3=sum(animalmodel{animalnum}.lmodel.rep{rep}.part{p}.nogo.X{lag}...
    (find(animalmodel{animalnum}.lmodel.rep{rep}.part{p}.nogo.Y{lag}==type),:),1);
lag=4;
b4=sum(animalmodel{animalnum}.lmodel.rep{rep}.part{p}.nogo.X{lag}...
    (find(animalmodel{animalnum}.lmodel.rep{rep}.part{p}.nogo.Y{lag}==type),:),1);
lag=5;
b5=sum(animalmodel{animalnum}.lmodel.rep{rep}.part{p}.nogo.X{lag}...
    (find(animalmodel{animalnum}.lmodel.rep{rep}.part{p}.nogo.Y{lag}==type),:),1);
lag=6;
b6=sum(animalmodel{animalnum}.lmodel.rep{rep}.part{p}.nogo.X{lag}...
    (find(animalmodel{animalnum}.lmodel.rep{rep}.part{p}.nogo.Y{lag}==type),:),1);
lag=7;
b7=sum(animalmodel{animalnum}.lmodel.rep{rep}.part{p}.nogo.X{lag}...
    (find(animalmodel{animalnum}.lmodel.rep{rep}.part{p}.nogo.Y{lag}==type),:),1);
lag=8;
b8=sum(animalmodel{animalnum}.lmodel.rep{rep}.part{p}.nogo.X{lag}...
    (find(animalmodel{animalnum}.lmodel.rep{rep}.part{p}.nogo.Y{lag}==type),:),1);

out=[out;sum([a1;a2;a3;a4;a5;a6;a7;a8],1);sum([b1;b2;b3;b4;b5;b6;b7;b8],1)];

end

inds = find(~isnan(sum(out,1)));