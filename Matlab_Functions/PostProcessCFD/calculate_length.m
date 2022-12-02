function [x_net] = calculate_length(coord, idx)

% Extract axial coordinate
x = coord(:,1);

% Number of clusters
k = max(idx);

% Partition cell coordinates
x_clust = clustering(x, idx);

x_max = zeros(k,1);
for j = 1 : k
    x_max(j) = max(x_clust{j});
end

% Sort x max
x_net = sort(x_max, 'ascend');



    
    





end

