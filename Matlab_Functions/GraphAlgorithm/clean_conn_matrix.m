function [global_connect_new] = clean_conn_matrix(global_connect)
% global_connect = matrix with the connections of the cell on the first two
% columns. We need to eliminate duplicates, loops in the graph, boundary
% cells

% Remove boundaries
global_connect_no_bound = remove_boundary(global_connect);

% Remove any duplicate row
global_connect_new = unique(global_connect_no_bound, 'rows', 'stable');  





end

