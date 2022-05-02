% hazard rate
sample = [];
sample_cut = [];
sample_cut2 = [];
sample_g = [];
while length(sample)<1000
    
    st = exprnd(4);
    if 1% st<=10
        sample(end+1)  = st;
    end
end


while length(sample_cut)<1000
    
    st = exprnd(4);
    if  st<=10
        sample_cut(end+1)  = st;
    end
end

while length(sample_cut2)<1000
    
    st = exprnd(4);
    if  st<=15
        sample_cut2(end+1)  = st;
    end
end

while length(sample_g)<1000
    
      sample_g(end+1)  = normrnd(4,1);
end
%figure;histogram(sample,100,'Normalization','pdf')
[f,x] = ecdf(sample,'Function','cumulative hazard');
[fc,xc] = ecdf(sample_cut,'Function','cumulative hazard');
[fc2,xc2] = ecdf(sample_cut2,'Function','cumulative hazard');
[fcg,xcg] = ecdf(sample_g,'Function','cumulative hazard');
figure;plot(x,f);hold on; plot(xc,fc);hold on; plot(xc2,fc2);hold on; plot(xcg,fcg);



