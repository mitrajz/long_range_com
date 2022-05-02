function do_plot_vecs(vec_n_go,vec_n_nogo)
nlags = 8;
figure;
for i=1:nlags
    hold on; plot((vec_n_go{i}),'g');
    hold on; plot((vec_n_nogo{i}),'r');
end