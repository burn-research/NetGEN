function [C] = get_cluster_centroids(X, idx)
%
% function [C] = get_cluster_centroids(X, idx)
%
% This function will calculate the centroids given the data matrix X and
% the clustering index idx
% OUTPUT
%   C = k x nv vector of centroids

% Partition the data
X_clust = clustering(X, idx);

% Number of variables and number of clusters
k = length(X_clust);
[np, nv] = size(X);

% Calculate centroids
C = zeros(k, nv);
for i = 1 : k
    C(i,:) = mean(X_clust{i});
end



end

