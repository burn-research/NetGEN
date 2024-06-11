function [idx_new, k_new] = match_nr(cell_id, G, idx, nr, opt)
% This function reassign the clusters to match the user specified number of
% reactors. We can establish different criteria for the reassignment, such
% as we can start from the smallest and agglomerate it to its neighbour...
% we will see

% Initialize the new idx vector
idx_new = idx;

% Partition the computational cells
cell_clust = clustering(cell_id, idx_new);

% Check for empty clusters
for j = 1 : length(cell_clust)
    if isempty(cell_clust{j})
        mess = append('Cluster ', num2str(j), ' is empty, idx is updated');
        disp(mess);
        pause;
        
        % Update idx by decreasing the number of clusters
        for i = 1 : length(idx_new)
            if idx_new(i) > j
                idx_new(i) = idx_new(i) - 1;
            end
        end
    end
end

% Update the new number of clusters, re-group the cells
k = max(idx_new);
cell_clust = clustering(cell_id, idx_new);
            
% Initial check of the connectivity
[~, n_sub_graphs] = check_connectivity(G, cell_clust);

k_new = k;

% Check if we need agglomeration or splitting
if k > nr
    disp('Number of clusters greater than desired reactor. Agglomeration is needed');
elseif k == nr
    disp('Number of clusters equals number of reactors. No reassignment is needed');
else
    disp('Number of clusters lower than number of reactors. Splitting is needed');
end

%% Check criterion for clusters re-assignment
if isfield(opt, 'ClusterReassignmentCriterion') == false
    warning('No criterion for cluster reassignment was specified. The volume criterion will be used by default');
    crit = 'volumes';
else
    crit = opt.ClusterReassignmentCriterion;
end

%% Start re-assigning reactors

if strcmp(crit, 'volumes') == true
    if isfield(opt, 'CellsVolume') == false
        error('Reassignment based on volume selected but CellsVolume is absent as field in opt');
    end

    disp('Reassigning clusters to match desired number of reactors...');

    while k_new > nr 

        % Calculate the volume of the clusters
        V_cells = opt.CellsVolume;                  % Vector of the volume of all cells
        V_clust = clustering(V_cells, idx_new);     % Volume of cells clustered in cell array
        Vi = zeros(k_new,1);
        for i = 1 : k_new
            Vi = sum(V_clust{i});
        end
        
        % Find smallest cluster
        [~,sid] = min(Vi);
    
        % Cluster cells
        cell_clust = clustering(cell_id, idx_new);
    
        % Scan through the cells in the subgraphs and find their
        % neighbors. Store the neighbors that do not belong to this
        % cluster in another vector
        neighb_clust = [];
        count = 0;
        cell_i = cell_clust{sid};
        for l = 1 : length(cell_i)
            neighb = neighbors(G, cell_i(l));           % Neighbors of each cells
            clust_neighb = idx_new(neighb);             % Neighbors clusters
    
            % Find external clusters
            ext_clust = clust_neighb(clust_neighb ~= sid);
    
            % If it is empty, keep idx_new untouched and proceed to the
            % next cell
            if isempty(ext_clust) == false
                neighb_clust(count+1:count+length(ext_clust)) = ext_clust;
                count = count + length(ext_clust);
            end
    
        end
    
        % Now re-assign the cluster based on the neighbor cluster
        % with more communicating cells
        neighb_clust = neighb_clust';
        if isempty(neighb_clust) == true
            warning('A cluster seems to not have any neighbour. Please check.');
        end
    
        % Count the neighbor clusters
        [clust_count, clust_index]  = groupcounts(neighb_clust);
    
        % Get the neighbour cluster with more neighbour cells
        [~,idmax] = max(clust_count);
        kmax = clust_index(idmax);
    
        % Update
        idx_new(idx_new==sid) = kmax;
        idx_new(idx_new>sid)  = idx_new(idx_new>sid)-1;
        k_new = max(idx_new);

    end

    disp('Reassignment Done');

end








end

