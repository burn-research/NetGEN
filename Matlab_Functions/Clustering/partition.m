function [nz_X_k, nz_idx_clust, k] = partition(X, idx, k, min_var)
%
% function [nz_X_k, nz_idx_clust, k] = partition(X, idx, k, min_var)
%
% This subroutine performs the partition the data into clusters


if nargin == 3
    [rows columns] = size(X);

    idx_clust = cell(k, 1);
    n_points = zeros(k, 1);    
    for j = 1 : k
        idx_clust{j} = find(idx == j);
        n_points(j) = size(idx_clust{j}, 1);
        if (n_points(j) < columns/10)
            fprintf('\nNo points in the cluster n. %d, cluster removed \n', j);
        end
    end
    nz_idx = find(n_points > columns);
    k_new = size(nz_idx, 1);
    k = k_new;
    nz_X_k = cell(k, 1);
    nz_idx_clust = cell(k, 1);
    for j = 1 : k
        nz_X_k{j} = zeros(n_points(j), columns);
        nz_idx_clust{j} = idx_clust{nz_idx(j)};
        nz_X_k{j} = X(nz_idx_clust{j}, :);
    end
    
else
    [rows columns] = size(X);

    idx_clust = cell(k, 1);
    n_points = zeros(k, 1);
    for j = 1 : k
        idx_clust{j} = find(idx == j);
        n_points(j) = size(idx_clust{j}, 1);
        if (n_points(j) < columns/10 || min_var(j) < 0.9)
            fprintf('\nNo points in the cluster or variance too low n. %d, cluster removed \n', j);
        end
    end
    
    nz_idx = find(n_points > columns);
    nz_idx_var = find(min_var > 0.01);
    
    % Take the one with minimum length
    if length(nz_idx) < length(nz_idx_var)
        k_new = size(nz_idx, 1);
        k = k_new;
        nz_X_k = cell(k, 1);
        nz_idx_clust = cell(k, 1);
        for j = 1 : k
            nz_X_k{j} = zeros(n_points(j), columns);
            nz_idx_clust{j} = idx_clust{nz_idx(j)};
            nz_X_k{j} = X(nz_idx_clust{j}, :);
        end
        
    else
        k_new = size(nz_idx_var, 1);
        k = k_new;
        nz_X_k = cell(k, 1);
        nz_idx_clust = cell(k, 1);
        for j = 1 : k
            nz_X_k{j} = zeros(n_points(j), columns);
            nz_idx_clust{j} = idx_clust{nz_idx_var(j)};
            nz_X_k{j} = X(nz_idx_clust{j}, :);
        end
        
    end
        
    
end
    


 
       