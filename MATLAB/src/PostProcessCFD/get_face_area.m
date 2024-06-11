function [A] = get_face_area(data_boundaries, idx)
% This function creates a vector containing the boundary area of the cells
% in m2
% INPUTS:
%   data_boundaries = data from Boundary_cells (output file from fluent)
%   idx = clustering label vector
%
% OUTPUTS: 
%   A = vector with the face area of each reactor

nc = length(idx);               % Number of cells
Af = zeros(nc,1);               % Initialize vector of face areas
bc = data_boundaries(:,1);      % Boundary cell ids
ba = data_boundaries(:,4);      % Boundary cell face area

for i = 1 : nc 
    % Search if cell i is boundary cell
    id = find(bc == i);
    if isempty(id) == true
        Af(i) = 0;
    else
        Af(i) = ba(id(1));
    end
end

% Get face area of the reactors
A_clust = clustering(Af, idx);  % Cell type array with cell areas in each cluster
nr = max(idx);                  % Number of reactors
A = zeros(nr,1);
for i = 1 : nr
    A(i) = sum(A_clust{i});
end





end