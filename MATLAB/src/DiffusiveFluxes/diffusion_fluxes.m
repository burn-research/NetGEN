function [Dm, Ab, DD] = diffusion_fluxes(data_all, G, idx, opt)
%
% function [Dm] = diffusion_fluxes(data_all, G, idx)
%
% This function is used to calculate the estimated diffusive fluxes across
% the various clusters (or reactors) via simplified assumption:
%
%   INPUTS:
%
%       data_all = struct file with all the data imported (must have
%       Viscosity and Connectivity and Solution)
%
%       G = graph object representing the computational mesh (See
%       create_graph.m)
%
%       idx = clustering labelling vector
%
%       opt = opt struct file with options
%
%   OUTPUTS:
%
%       Dm = matrix of diffusive fluxes (size Nr*Nr). The element ij-th is
%       the diffusive flux from reactor i to reactor j in kg/s
%
%       Ab = boundary areas between the clusters (m2)
%
%       DD = geometric distance between the clusters (m)

% Extract data
coord = data_all.Coordinates;
Y     = data_all.Solution(:,2:end); % Take only species

DD = calc_cluster_distance(coord, idx);

% Check dimensions
if length(Y) ~= length(idx)
    error('Check dimensions of inputs. Dimension of Y and idx must agree');
end

% Vector of nodes
nodes = [1:1:length(idx)]';

% Subgraphs
[~, ~, ~, ~, ~, subgraphs] = ...
    connectivity_study(G, nodes, idx, false, coord);

% Calculate contact surface
val_neighbors = data_all.Connectivity;
Ab = calc_boundary_surface(G, subgraphs, idx, nodes, val_neighbors);
Ab = Ab * 2 * pi;

% Calculate diffusion fluxes
if isfield(data_all, 'Viscosity') == false
    error('Effective viscosity must be provided');
end

% Check if mass or mole fraction are provided
if isfield(opt, 'Basis') == false
    warning('The basis (mol or mass) was not provided. Mass fractions will be assumed by default');
else
    basis = opt.Basis;
    if strcmp(basis, 'mol') || strcmp(basis, 'mole') || strcmp(basis, 'molar')
        sp_labels = data_all.SpeciesLabels;
        sp_labels_c = rewrite_fluent_labels(sp_labels);
        Y = mole_to_mass(Y, sp_labels_c);
        fprintf('Mole basis was specified so species were converted to mass fractions');
    end
end

Dm = estimate_diffusion_fluxes(Y, DD, Ab, data_all.Viscosity, idx);
disp('Diffusive fluxes successfully estimated');




end

