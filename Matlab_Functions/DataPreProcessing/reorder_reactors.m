function [idx_new] = reorder_reactors(idx, OrdVar)
% This function reorder the reactors according to the variable OrdVar. 

% Get number of clusters
k = max(idx);

% Cluster data and get mean variable in clusters
VClust_cell = clustering(OrdVar, idx);
VClust_arr = zeros(k,1);
for i = 1 : k
    VClust_arr(i) = max(VClust_cell{i});
end

% Sort array
[~, id_sort] = sort(VClust_arr, 'ascend');

% Update idx_new
idx_new = idx;
for i = 1 : k
    idx_new(idx==id_sort(i)) = i;
end











end

