function [f] = get_mixture_fraction(data_solution, fuel)
% This function returns the bilger mixture fraction given the data
% structure field data_solution and the fuel structure containing the
% composition of the fuel. The data structure field must contain the
% chemical species and their labels. The chemical species must be given as
% mole fractions

% Extract data and labels
data = data_solution.data;
labels = data_solution.textdata;

% Check if is a 2D or 3D case
if data(1,4) > 100
    dim = 2;
else
    dim = 3;
end

% Get the corresponding data and labels
data = data(:,dim+2:end);
labels = labels(dim+2:end);

% Get mixture fraction
sp_labels_c = rewrite_fluent_labels(labels);
Y = mole_to_mass(data, sp_labels_c);
f = mixture_fraction(Y, labels, fuel);





end

