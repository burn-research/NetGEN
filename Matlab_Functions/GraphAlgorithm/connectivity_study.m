function [n_subgraphs, subgraph_sizes, nodes_in_subgraphs, spared_cells, spared_cells_tot, subgraphs] = connectivity_study(G, cell_id, idx, plot, coord)
% This function will assess useful parameters to evaluate the geometrical
% connectivity of a clustering. 
%   INPUTS:
%       G = graph of the mesh (see the function create_graph.m)
%       cell_id = vector of cell id number
%       idx = clusters labels vector
%
%   OUTPUTS:
%       n_sub_graphs = vector of the total number of spared graphs in each
%                       cluster
%
%       sub_graph_sizes = cell array with k vector of the sizes of each
%       spared graph in each cluster
%
%       nodes_in_subgraphs = cell array with k elements representing the
%       partition of the nodes in each subgraph, for each cluster
%
%       spared_cells = vector with spared cells in each cluster


%% Partition the cell and the coordinates
cell_partition = clustering(cell_id, idx);
coord_partition = clustering(coord, idx);

%% Evaluate cluster sizes
k = max(idx);

clust_size = zeros(k,1);
for j = 1 : k
    clust_size(j) = length(cell_partition{j});
end

%% Perform connectivity analysis
nodes_in_subgraphs = cell(k, 1);        % Cell array with nodes in each subgraphs, for each cluster
subgraph_sizes = cell(k, 1);            % Cell array with subgraphs sizes, for each cluster
spared_cells = zeros(k,1);              % Number of spared cell in each cluster
n_subgraphs = zeros(k,1);               % Number of subgraphs in each cluster
subgraphs = cell(k, 1);
for i = 1 : k
    H = subgraph(G, cell_partition{i}) ;                    % Subgraph of cluster j
    [bins, binsizes] = conncomp(H, 'OutputForm', 'cell');   % Bins is a cell array that contains nodes in each subgraphs
    n_subgraphs(i) = length(binsizes);                      % Number of subgraphs in cluster j
    nodes_in_subgraphs{i} = bins;                           % Cell array with node id in each subgraph
    subgraph_sizes{i} = binsizes;                           % Size of each subgraph
    binsize_sort = sort(binsizes, 'descend');
    spared_cells(i) = sum(binsize_sort(2:end));             % Total number of spared cell in cluster j
    subgraphs{i} = H;
end

spared_cells_tot = sum(spared_cells);                      % Total number of spared cells
mess = sprintf('Total number of spared cells = %d', spared_cells_tot);
disp(mess);


%% Plot the results
if plot == true
    subplot(1,2,1);
    bar(n_subgraphs);
    xlabel('Number of cluster');
    ylabel('Number of subgraphs');
    title('Number of subgraphs in each cluster');
    
    figure; subplot(1,2,2);
    bar(spared_cells);
    xlabel('Number of cluster');
    ylabel('Number of spared cells');
    title('Number of spared cells in each cluster');
end

