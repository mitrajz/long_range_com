
v1n = [];
lmn = [];
load('cells_FF_pl20_an150_lw20_exG1_onlyC1_onlyS0_plstyle1.mat')
for i=unique(cellfun(@(x) x.simulcode,V1cells))
    v1n=[v1n; numel(find(cellfun(@(x) x.simulcode,V1cells) ==i))];
    lmn=[lmn; numel(find(cellfun(@(x) x.simulcode,LMcells) ==i))];
end

load('cells_FB_pl20_an150_lw20_exG1_onlyC1_onlyS0_plstyle1.mat')
for i=unique(cellfun(@(x) x.simulcode,V1cells))
    v1n=[v1n; numel(find(cellfun(@(x) x.simulcode,V1cells) ==i))];
    lmn=[lmn; numel(find(cellfun(@(x) x.simulcode,LMcells) ==i))];
end

% rounded numbers:
[round(mean(v1n)),(round(std(v1n)))] % % mean across animals, std across animals
[round(mean(lmn)),(round(std(lmn)))] % % mean across animals, std across animals

sum(v1n)
sum(lmn)

sum(v1n)+sum(lmn)
%%
uv1n = [];
ulmn = [];
load('cells_U_FB_pl20_an150_lw20_exG0_onlyC0_onlyS0_plstyle1.mat')
for i=unique(cellfun(@(x) x.simulcode,V1cells))
    uv1n=[uv1n; numel(find(cellfun(@(x) x.simulcode,V1cells) ==i))];
    ulmn=[ulmn; numel(find(cellfun(@(x) x.simulcode,LMcells) ==i))];
end

sum(uv1n)
sum(ulmn)
sum(uv1n)+sum(ulmn)

