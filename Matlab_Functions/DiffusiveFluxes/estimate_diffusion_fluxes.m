function [Dm] = estimate_diffusion_fluxes(Y, Dc, Ab, viscosity, idx)
% This function estimates the diffusive fluxes from the clustered solution
% for each species, from each reactor to the neighboring reactors.
% INPUTS:
%       Y     = mass fraction matrix (np x ns)
%       Dc    = distance between clusters (k x k)
%       Ab    = surface of contact between clusters (k x k)
%       viscosity = effective viscosity of cells (np x 1)
%       idx = clustering labelling vector (integers np x 1)
%
% OUTPUTS:
%       D = cell array of diffusive fluxes ns x (k x k)

%% Dimensions checks and preliminary variables
[np, ns] = size(Y);

if np ~= length(viscosity) || np ~= length(idx)
    error('Please check dimensions in input');
end

% Number of clusters
k = max(idx);

if k ~= size(Dc, 1) || k ~= size(Ab, 1)
    error('Number of cluster must agree with size of Dc and Ab');
end

%% Calculate viscosity of the clusters
visc_clust = clustering(viscosity, idx);
mu_c = zeros(k, 1);
for i = 1 : k
    mu_c(i) = mean(visc_clust{i});
end

%% Calculate diffusive fluxes
Dm = cell(ns, 1);
for l = 1 : ns
    Yclust = clustering(Y(:,l), idx);    
    Di = zeros(k, k);
    for i = 1 : k
        for j = 1 : k

            if i ~= j
                Di(i,j) = Ab(i,j) * (mu_c(j)/0.7) * (mean(Yclust{j}) - mean(Yclust{i}))/Dc(i,j);
            end
        end
    end

    Dm{l} = Di;
end











end

