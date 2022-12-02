function [idx_new, k_new] = reassign_nodes(cell_id, G, idx, opt)
% This function will re-label the computational cells labeled by the
% clustering vector idx, embedded in the graph G, which contains the
% information of the connectivity. The cell will be re-assigned based on
% subgraphs similarities in order to ensure the connectivity of the domain
% INPUTS
%   cell_id = vector containing cell identification number (from Fluent)
%       with dimension n_cells x 1
%   
%   G = graph representing the mesh (see create_graph.m)
%
%   idx = clustering label vector with dimension n_cells x 1
%   
%   coord = matrix containing cell coordinates with dimensions n_cells x
%   dim (dim is the geometric dimension of the problem)
%
% OUTPUT
%
%   idx_new = new clustering label vector
%
%   k_new = new number of clusters

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

% Set parameters for routine
itmax = 1000;
it = 1;
counter = 0;

% Parameter initialization
convergence = false; 
k_new = k;

while it < itmax && convergence == false
    
    % Scan through the clusters
    for j = 1 : k_new
        cell_clust_j = cell_clust{j};   % Cell ids contained in cluster j
        
        % Check if is an empty cluster
        if isempty(cell_clust_j)
            mess = append('Cluster ', num2str(j), ' is empty, it will be deleted');
            disp(mess);
            
            % If the cluster is empty, then the idx vector is updated by
            % subtracting 1 to all the clusters with a label bigger than j
            for i = 1 : length(idx_new)
                if idx_new(i) > j
                    idx_new(i) = idx_new(i) - 1;
                end
            end
            
            % Update the number of cluster
            k_new = k - 1;
            break
        
        % If the cluster is not empty
        else
        
            % Build the subgraph
            H = subgraph(G, cell_clust_j);

            % Evaluate the disconnected components
            % This is a cell type array with dimension n_subgraph x 1. Each
            % element contain the ids of the cell in the subgraph
            bins = conncomp(H, 'OutputForm', 'cell');                                    
            binsize = cellfun('length', bins);  % Subgraph dimensions (number of cells in each subgraph)

            % Sort the disconnected graphs according to their size, from the
            % biggest to the smallest
            [binsize, ind] = sort(binsize, 'descend');
            bins = bins(ind);

            % Now look for small subgraphs
            if length(binsize) > 1

                reassign_graph = false;
                if isfield(opt, 'ReassignCriterion') == false
                    crit = 'n_cells';
                else
                    crit = opt.ReassignCriterion;
                end

                % Criterion based on n_cells
                if strcmp(crit, 'n_cells')
                    if binsize(2) < 1000
                        reassign_graph = true;
                    end

                % Calculate the volume of the small cluster found by the
                % cells in bins{2}
                elseif strcmp(crit, 'volume')
                    cell_i = zeros(binsize(2),1);
                    for l = 1 : binsize(2)
                        cell_i(l) = str2double(bins{2}{l});
                    end

                    if isfield(opt, 'CellsVolume') == false
                        error('Cells volume not specified');
                    end

                    V_cells = opt.CellsVolume;      % Vector of the volume of all cells
                    V_i = sum(V_cells(cell_i));     % Volume of the small cluster

                    if isfield(opt, 'VolumeThreshold') == false
                        fprintf('Volume threshold not specified, 0.01 will be used as default \n');
                        opt.VolumeThreshold = 0.01;
                    end

                    if V_i < opt.VolumeThreshold*sum(V_cells)
                        reassign_graph = true;
                        fprintf('Small cluster in volume');
                    end
                end

                
                % If the graph is smaller than 100 cells
                if reassign_graph == true
                    disp('Small graph found, it will be reassigned');
                    cell_i = zeros(binsize(2), 1);  
                    for l = 1 : binsize(2)
                        cell_i(l) = str2double(bins{2}{l});
                    end

                    % Scan through the cells in the subgraphs and find their
                    % neighbors. Store the neighbors that do not belong to this
                    % cluster in another vector
                    neighb_clust = [];
                    count = 0;
                    for l = 1 : length(cell_i)
                        neighb = neighbors(G, cell_i(l));           % Neighbors of each cells
                        clust_neighb = idx_new(neighb);             % Neighbors clusters

                        % Find external clusters
                        ext_clust = clust_neighb(clust_neighb ~= j);

                        % If it is empty, keep idx_new untouched and proceed to the
                        % next cell
                        if isempty(ext_clust) == true
                            disp('No external clusters');

                        % If it is not empty store the cluster label in the vector
                        % and update the counter
                        else
                            neighb_clust(count+1:count+length(ext_clust)) = ext_clust;
                            count = count + length(ext_clust);
                        end

                    end

                    % Now re-assign the cluster based on the neighbor cluster
                    % with more communicating cells
                    neighb_clust = neighb_clust';

                    % Count the neighbor clusters
                    clust_count  = groupcounts(neighb_clust);

                    % Get the clusters that communicate with the subgraph
                    clust_unique = sort(unique(neighb_clust), 'ascend');

                    % Select criterion to re-assign the cells
                    if isfield(opt, 'Criterion') == false
                        opt.Criterion = 'Neighbor';
                    end

                    if strcmp(opt.Criterion, 'Neighbor') == true
                        % Re-assign the cells to the neighbor cluster with more
                        % neighboring cells
                        [~, id] = max(clust_count);
                         new_clust = clust_unique(id);
    
                        if isempty(new_clust)
                            warning('New cluster does not exist');
                        else
                            % Re-assign the cells and re-evaluate the connectivity of the
                            % clustering
                           idx_new(cell_i) = new_clust;
                        end
                        
                        counter = counter + 1;
                        mess = append(num2str(counter), ' subgraph moved');
                        disp(mess);
                        break
                    elseif strcmp(opt.Criterion, 'Euclidean') == true

                        if isfield(opt, 'X') == false
                            error('Please specify the X matrix via opt.X as input for the function');
                        else
                            X = opt.X;  % Data matrix
                            C = get_cluster_centroids(X, idx); % Cluster centroids

                            % Get centroids of the small subgraph
                            Ci = mean(X(cell_i));

                            % Get the distance between Ci and only the
                            % adjacent clusters
                            D = zeros(length(clust_unique),1);  % Array of distances
                            for s = 1 : length(clust_unique)
                                D(s) = sum((Ci - C(clust_unique(s))).^2);
                            end

                            % Find minimum distance
                            [~,id_min] = min(D);
                            new_clust = clust_unique(id_min);   % New cluster to re-assign the cell to
                            if isempty(new_clust)
                                warning('New cluster does not exist');
                            else
                                % Re-assign the cells and re-evaluate the connectivity of the
                                % clustering
                               idx_new(cell_i) = new_clust;
                            end
                            counter = counter + 1;
                            mess = append(num2str(counter), ' subgraph moved with Euclidean criterion');
                            disp(mess);
                            break
                        end
                    end

                % If it is a big cluster, split it   
                else
                    disp('Big cluster found, the original cluster will be splitted');
                    cell_i = zeros(binsize(2), 1);
                    k_new = k_new + 1;
                    for l = 1 : binsize(2)
                        cell_i(l) = str2double(bins{2}{l});
                    end
                    
                    
                    % Change the clustering index of the cells
                    for l = 1 : length(cell_i)
                        idx_new(cell_i(l)) = k_new;
                    end
                    
                    counter = counter + 1;
                    mess = append(num2str(counter), ' subgraph moved');
                    disp(mess);
                    break
                        

                end

            % If the cluster is composed by only one cell, move the cell and delete the cluster    
            elseif length(binsize) == 1 && binsize(1) == 1
                disp('Cluster with only one cell found, it will be reassigned');
                cell_i = str2double(bins{1}{1});
                neighb = neighbors(G, cell_i);
                clust_neighb = idx_new(neighb);
                new_clust = clust_neighb(1);
                
                % Update the cluster of the cell
                idx_new(cell_i) = new_clust;
                break
                
            end
                            
        end

    end
    
    % Re-evaluate the connectivity
    cell_clust_new = clustering(cell_id, idx_new); 
    [all_connected, n_sub_graphs_new] = check_connectivity(G, cell_clust_new);
    
    
    n_sub_graphs_new = sum(n_sub_graphs_new);
    n_sub_graphs = sum(n_sub_graphs);    

    if n_sub_graphs_new > n_sub_graphs
        warning('Number of disconnected graphs is raising!!!');
    elseif n_sub_graphs_new < n_sub_graphs
        disp('The connectivity has been improved');
    elseif n_sub_graphs_new == n_sub_graphs
        disp('Connectivity still the same');               
    end
    
    % Print actual number of clusters
    fprintf('\n Actual number of clusters: %d \n', k_new);
    
    % Check convergence
    if all_connected == true
        convergence = true;
        disp('The nodes has been reassigned and the clusters are now connected');
    else
        n_sub_graphs = n_sub_graphs_new;
        cell_clust = clustering(cell_id, idx_new);
        it = it + 1;
        if it == itmax
            disp('Maximum number of iterations reached and the domain is still not connected');
        end
    end
        
end

% % Ask to display the final results
% display = input('Enter yes to display the results, no to proceed without: ', 's');
% switch display
%     case 'yes'
%         H = plot_subgraphs(cell_id, idx_new, G, coord);
%     case 'no'
%         disp('Results not displayed');
% end        


end

