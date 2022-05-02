function y = convert_tableunit_to_templateunit(x,path)

[tableSpikeData, tableClusterData, strParams] = loadKiloSortedSpikes(path);
a=(sort(unique(tableSpikeData.spike_clusters)));


clusters=readNPY([path,'spike_clusters.npy']);
b=unique(clusters);

y = nan(size(x));
y(find(~isnan(x))) = a(x(find(~isnan(x))));
