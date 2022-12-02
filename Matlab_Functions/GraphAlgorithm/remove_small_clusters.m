function [idx_new, new_volumes] = remove_small_clusters(net_volumes, vol_thresh, idx, G)
% This function will remove clusters whose volume is below the threshold
% vol_thresh. It will produce a new cluster labelling vector idx_new by
% updating the current idx. The clusters will be reassigned according to
% the neighbor with more communicating cells

% Number of clusters
k = max(idx);
k_new = k;

% Number of cells
nc = length(idx);
cell_id = [1:1:nc]';

% Partition the cells
cell_part = clustering(cell_id, idx);

% New volumes vector (new volumes should be a column vector)
[n,m] = size(net_volumes);
if n == 1 && m > 1
    new_volumes = net_volumes';
elseif m == 1 && n > 1
    new_volumes = net_volumes;
else
    error('Check dimensions of net_volumes array, something is wrong');
end

% Initialize the new index vector
idx_new = idx;

while isempty(find(new_volumes < vol_thresh)) == false
    
    clust_to_remove = find(new_volumes < vol_thresh);
    
    for i = 1 : length(clust_to_remove)
        
        cell_clust = cell_part{clust_to_remove(i)};
        ext_neighb = [];
        count = 1;
        
        % Detect neighbors of those cells keep only the external neighbors
        % (different clusters)
        for j = 1 : length(cell_clust)
            neighb = neighbors(G, cell_clust(j));
            for l = 1 : length(neighb)
                if idx_new(neighb(l)) ~= clust_to_remove(i)
                    ext_neighb(count) = idx_new(neighb(l));
                    count = count + 1;
                end
            end
        end
        
        % Now count the elements
        [nn, nid] = groupcounts(ext_neighb');
        [~, max_id] = max(nn);
        new_clust = nid(max_id);                % New cluster
        
        % Print the message
        fprintf('Cluster %d has been moved to cluster %d \n', clust_to_remove(i), new_clust);
        
        % Assign the cells to the new cluster
        idx_new(cell_clust) = new_clust;
        
        % Calculate the new volumes
        new_volumes(new_clust) = new_volumes(new_clust) + new_volumes(clust_to_remove(i));
        
        % Set to zero the old volume
        new_volumes(clust_to_remove(i)) = 0;
        
        % Since the update to a cluster can modify the others, it is
        % necessary to exit from the for cycle and do one cluster at a time
        break
        
    end
    
    % Update the volumes by removing the zero elements
    vol_update = removerows(new_volumes, 'ind', clust_to_remove(1));
    new_volumes = vol_update;
    
    % Update the number of clusters
    k_new = k_new - 1;
    
    % Modify the idx of the clusters above the one deleted
    for j = 1 : length(idx)
        if idx_new(j) > clust_to_remove(1)
            idx_new(j) = idx_new(j) - 1;
        end
    end
    
    % Update cell partition
    cell_part = clustering(cell_id, idx_new);
    
    
end

% Check connection again
cell_clust = clustering(cell_id, idx_new);
[all_connected, n_sub_graph] = check_connectivity(G, cell_clust);

if all_connected == false
    warning('Clusters are not connected');
end

end


    
    
    
       
        
        
        
        
        
        

