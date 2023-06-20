close all
clear all

set(0, 'defaulttextfontsize', 12);
set(0, 'defaultaxesfontsize', 12);
set(0, 'defaulttextfontname', 'Times');
set(0, 'defaultfigureunits', 'centimeters');

%% Test distance
% Import data
data = importdata('/Users/matteosavarese/Desktop/Dottorato/Data_Furnace/Data_CFD/2D/25mm/50CH4_50H2_phi08_kee/data_solution');
val = data.data;
coord = val(:,2:3);
Y = val(:,5:end);

idx = kmeans(coord, 10);

DD = calc_cluster_distance(coord, idx);
scatter(coord(:,1), coord(:,2), 1, idx, 'filled');

%% Connectivity study

% Import mass flow data
data_mass = importdata('/Users/matteosavarese/Desktop/Dottorato/Data_Furnace/Data_CFD/2D/25mm/50CH4_50H2_phi08_kee/Neighbours_cell_flows');

% Remove the boundary cells from the connectivity matrix (refers to remove_boundary.m)
data_mass = data_mass.data;
data_mass = data_mass(:,4:end);
mass_data = remove_boundary(data_mass);

% Build the graph (refers to create_graph.m)
nodes = [1:1:length(idx)]';
conn_matrix = mass_data(:,1:2);
G = create_graph(nodes, conn_matrix);

% Analyze the connectivity (refers to check_connectivity.m and clustering.m)
% Partition the cell id's
cell_clust = clustering(nodes, idx);

% Check for empty clusters
idx_new = idx;
for j = 1 : length(cell_clust)
    if isempty(cell_clust{j})
        mess = append('Cluster ', num2str(j), ' is empty');
        disp(mess);
        for i = 1 : length(idx)
            if idx_new(i) > j
                idx_new(i) = idx_new(i) - 1;
            end
        end
    end
end

% Evaluate the new number of clusters
k_new = max(idx_new);
idx = idx_new;
mess = append('Now we have ', num2str(k_new), ' clusters');
disp(mess);

% Re-group the cells
cell_clust = clustering(nodes, idx_new);

% Check connectivity
[all_connected, n_sub_graph] = check_connectivity(G, cell_clust);

% Get graph and subgraphs
plt = false;
[n_subgraphs, subgraph_sizes, nodes_in_subgraphs, spared_cells, spared_cells_tot, subgraphs] = ...
    connectivity_study(G, nodes, idx, plt, coord);

% Calculate area
data_volumes = importdata('/Users/matteosavarese/Desktop/Dottorato/Data_Furnace/Data_CFD/2D/25mm/50CH4_50H2_phi08_kee/data_volumes');
areas = data_volumes.data(:,4);

% Neighbor data
data_neighbors = importdata('/Users/matteosavarese/Desktop/Dottorato/Data_Furnace/Data_CFD/2D/25mm/0CH4_100H2_phi08_kee/Neighbours_cell_flows');
val_neighbors = data_neighbors.data;

Ab = calc_boundary_surface(G, subgraphs, idx, nodes, val_neighbors);

sc = scatter(coord(:,1), coord(:,2), 10, idx, 'filled');
fig = gcf; fig.Position = [15 15 16 8];

%% Estimate diffusion fluxes

data_visc = importdata('/Users/matteosavarese/Desktop/Dottorato/Data_Furnace/Data_CFD/2D/25mm/50CH4_50H2_phi08_kee/data_viscosity');
viscosity = data_visc.data(:,4);
[Dm] = estimate_diffusion_fluxes(Y, DD, Ab, viscosity, idx);

%% Test neighbor searches
data_neighbors = importdata('/Users/matteosavarese/Desktop/Dottorato/Data_Furnace/Data_CFD/2D/25mm/0CH4_100H2_phi08_kee/Neighbours_cell_flows');
val_neighbors = data_neighbors.data;

id1 = val_neighbors(:,4);
id2 = val_neighbors(:,5);

