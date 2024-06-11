% Condition.m - This function perform data conditioning
%
% function [bin_data, idx_clust] = condition(data, z, n_bins, min_z, max_z,
% z_stoich) 
%
% INPUTS
% 
% data            = Matrix of state variables
% z               = Conditioning variable vector
% n_bins          = Number of bins to condition the data into
% OPTIONAL:
% min_z           = Minimum value of the coinditioning variable (if not
%                   selected the minimum absolute value of z will be used)
% max_z           = Maximum value of the coinditioning variable (if not
%                   selected the maximum absolute value of z will be used)
% z_stoich        = If conditioning on mixture fraction, the value of the
%                   stoichiometric mixture fraction can be provided. If
%                   not, a uniform partition (into uniform bins) is
%                   performed   
% OUTPUTS
%
% bin_data        = Cell matrix od conditioned data. Each cell stores the
%                   data corresponding to a single bin
% idx_clust       = Cell vector storing the indexes of the original data
%                   points in the bins. 

function [bin_data, idx_clust] = condition(data, z, n_bins, min_z, max_z, ...
    z_stoich)

% Check the inputs

if size(z, 2) > 1
    error('The conditioning variable MUST be a vector');
end

if size(data, 1) ~= size(z, 1)
    error('Inputs do not have the same length')
end

if nargin == 3
    min_z = input('Enter the min value for the conditioning variable: ');
    if (isempty(min_z))
        fprintf( '\nYou did not specify min_z. The min value of z will be selected \n');
        min_z = min(z);
    end
    max_z = input('Enter the max value for the conditioning variable: ');
    if (isempty(max_z))
        fprintf( '\nYou did not specify max_z. The max value of z will be selected \n');
        max_z = max(z);
    end
    z_stoich = input('If consitioning on mixture fraction, enter the stoichiometric mixture fraction ');
    if (isempty(z_stoich))
        fprintf( '\nYou did not specify z_stoich. Data will be conditioned into uniform bins in z space \n');
        cond_type = 'uniform';
    else
    cond_type = 'non_uniform';
    end
end

if nargin == 4
    error('You must specify a val;ue for z_max');
end

if nargin == 5
    z_stoich = input('If consitioning on mixture fraction, enter the stoichiometric mixture fraction ');
    if (isempty(z_stoich))
        fprintf( '\nYou did not specify z_stoich. Data will be conditioned into uniform bins in z space \n');
        cond_type = 'uniform';
    else
    cond_type = 'non_uniform';
    end
end

if nargin == 6
    cond_type = 'non_uniform';
end

% Number of intervals

n = n_bins + 1;

% Partition of the space

if n_bins == 1
    ints = linspace(min_z, max_z, n);
elseif strcmp(cond_type, 'uniform')
    ints = linspace(min_z, max_z, n);
elseif strcmp(cond_type, 'non_uniform')
    ints_1 = linspace(min_z, z_stoich, ceil(n/2));
    ints_2 = linspace(z_stoich, max_z, ceil((n+1)/2));
    ints = [ints_1(1:ceil(n/2-1)) ints_2];
end

% Bin data matrices initialization

bin_data = cell(n_bins, 1);
idx_clust = cell(n_bins, 1);

%  Partition

for bin = 1 : n_bins  
    idx_clust{bin} = find((z >= ints(bin)) & (z <= ints(bin+1)));
    bin_data{bin} = data(idx_clust{bin}, :);
end