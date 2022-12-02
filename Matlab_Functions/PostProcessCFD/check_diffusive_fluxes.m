function [mf_diffusion] = check_diffusive_fluxes(diffusive_flux_file, idx)
% This function calculates the diffusive fluxes between clusters given the
% input file taken from Fluent. The output is a cell of size n_sp containing
% k x k martices with the mass diffusive fluxes of all the chemical species
    
%% Import the data
diff = importdata(diffusive_flux_file);
data_diff = diff.data;
data_diff_nobound = remove_boundary(data_diff);

% Number of species
[~,n_sp] = size(data_diff_nobound);
n_sp = n_sp-2;
mf_diffusion = cell(n_sp,1);

%% Calculate diffusive fluxes
% Initialize mass flowrate matrix. Mass flowrate (i,j) is the mass flowing
% from reactor i to j
for s = 1 : n_sp
    mf = zeros(k,k);
    mass_data = [data_diff_nobound(:,1:2) data_diff_nobound(:,s+2)];
    
    % Scanning through neighbor cells
    for j = 1 : size(mass_data, 1)
        id_1 = mass_data(j,1)+1;
        id_2 = mass_data(j,2)+1;
        
        % If neighbor cells belong to different cluster, compute the mass
        % flowrate
        if id_2 ~= 0 && id_2 ~= -1 && isnan(id_1) == false
            clust_1 = idx(id_1);
            clust_2 = idx(id_2);
            if clust_1 ~= clust_2
                mi = mass_data(j,3);
                
                % If mass flowrate > 0, the flow is going from clust_1 to
                % clust_2, so store it in mf(clust_1, clust_2); if
                % the mass flowrate < 0 is the opposite, store it in
                % mf(clust_2, clust_1);
                
                if mi > 0
                    mf(clust_1, clust_2) = mf(clust_1, clust_2) + mi;
                else
                    mf(clust_2, clust_1) = mf(clust_2, clust_1) - mi;
                end
            end
        end
    end
    mf_diffusion{s} = mf;
end
