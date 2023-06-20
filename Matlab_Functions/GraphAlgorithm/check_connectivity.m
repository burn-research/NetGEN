function [all_connected, n_sub_graphs] = check_connectivity(G, cell_partition)
%
% [all_connected, n_sub_graphs] = check_connectivity(G, cell_partition)
%
% This function will check if a clustering results in geometrically
% connected clusters.
%
% INPUTS:
%   G = graph object representing the CFD mesh
%
%   cell_partition = a cell type array containing the cell id's partition
%
% OUTPUTS:
%   all_connected = a true or false variable. True if all clusters are
%                   connected, false if at least one cluster is disconnected
%
%   n_sub_graphs = vector containing the number of disconnected components
%                  in each cluster

k = length(cell_partition);
n_sub_graphs = zeros(k,1);

% For each partition, extract the subgraph H from the main graph G and
% check how many disconnected components H has
for j = 1 : k
    cell_clust = cell_partition{j};         % Cell ids contained in cluster j
    H = subgraph(G, cell_clust);            % Subgraph representing cluster j
    [~, binsizes] = conncomp(H);            % Disconnected component in graph H
    n_sub_graphs(j) = length(binsizes);     % Number of disconnected components
end

all_connected = true;
for j = 1 : k
    if n_sub_graphs(j) ~= 1
        all_connected = false;
        mess = sprintf('Cluster %d is disconnected in %d components', j, n_sub_graphs(j));
        disp(mess);
    end
end

if all_connected == true
    disp('All the clusters are connected');
end

        
        
    

end

