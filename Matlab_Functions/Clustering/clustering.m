function [X_clust, index_clust] = clustering(X, idx)
% Partition the data matrix X in clusters according to the indexes idx

n_clust = max(idx);
X_clust = cell(n_clust,1);
index_clust = cell(n_clust,1);

for i = 1 : n_clust
    points = find(idx == i);
    X_clust{i} = X(points, :);
    index_clust{i} = points;
end



end

