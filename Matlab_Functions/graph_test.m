close all;
clear all;

% This script will test the function re-assign nodes on test graphs
cell_id = [1:1:20]';
idx = randi(5,20,1);

% Create the connectivity Matrix
s = [1 2 4 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20]';
t = [2 3 5 1 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 1]';

conn_matrix = [s t];

coord = randn(20,3);

G = create_graph(cell_id, conn_matrix);
plot(G);
plot_subgraphs(cell_id, idx, G, coord);
pause

% Re-assign nodes
idx_new = reassign_nodes(cell_id, G, idx, coord);

