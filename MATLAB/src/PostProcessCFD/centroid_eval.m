function [C_geom] = centroid_eval(coord, idx)
% This function calculate the physical centroids of a partition 
% vector idx 

[~, n_dim] = size(coord);
n_clust = max(idx);
coord_clust = clustering(coord, idx);

C_geom = zeros(n_clust, n_dim);
for j = 1 : n_clust
    if isempty(coord_clust{j}) == true
        C_geom(j,:) = zeros(1, n_dim);
    else
        C_geom(j,:) = mean(coord_clust{j});
    end
end

end

