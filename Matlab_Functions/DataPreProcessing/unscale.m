% This function provides a tool to rescale the data after PCA has been
% performed
%
% INPUTS:
% 
% scaled_data   = data matrix to be unscaled. Each row is an observation,
%                 each column a variable
% gamma         = scaling parameter adopted for each avriable (column)
%
%
% OUTPUT
%
% unscaled_data = Matrix of unscaled data

function [unscaled_data] = unscale(scaled_data, gamma)

[rows, columns] = size(scaled_data);

% Initialization

unscaled_data = zeros(rows, columns);

for j = 1 : columns
    unscaled_data(:, j) = scaled_data(:, j) .* gamma(j);
end

%save recovered_original_data unscaled_data