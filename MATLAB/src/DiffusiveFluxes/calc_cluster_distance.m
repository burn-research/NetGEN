function [DD] = calc_cluster_distance(coord, idx)
% This function returns the k x k martrix of characteristic distances
% between clusters (in the geometric coordinate space) given the
% coordinates vector coord and the clustering label vector idx

% Check shapes correctness
[np, ndim] = size(coord);
if np ~= length(idx)
    error('Length of coord and idx must agree');
end

% Get number of clusters
k = max(idx);

% Get the cluster centroids
Cgeom = zeros(k, ndim);
coord_clust = clustering(coord, idx);
for i = 1 : k
    Cgeom(i,:) = mean(coord_clust{i});
end

% Get the distance matrix
DD = zeros(k,k);
for i = 1 : k
    for j = 1 : k
        DD(i,j) = sum((Cgeom(i,:) - Cgeom(j,:)).^2)^0.5;
    end
end







end

