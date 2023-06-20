function [G] = create_graph(cell_id, conn_matrix)
%
% function [G] = create_graph(cell_id, conn_matrix)
%
% This function will create an output graph G from the vector of the cell
% node identification number and the connectivity matrix from Fluent

n = length(cell_id);   % Number of cells

% You have to give names to the cell nodes instead of simple numbers,
% otherwise every time a subgraph is extracted, it will be re-indexed from
% 1
names = cell(n,1);
for j = 1 : n
    names{j} = num2str(cell_id(j));   % Name each node as its cell id
end

% Clean the connectivity matrix removing cell from boundary faces (with no
% neighbor cell)
conn_matrix_no_bound = remove_boundary(conn_matrix);

% Create the graph
G = graph(conn_matrix_no_bound(:,1)', conn_matrix_no_bound(:,2)', [], names);

end

