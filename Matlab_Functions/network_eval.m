function [mass_corrected, net_volumes, bc_mass_flowrates, idx, X, infos] = network_eval(k, inlet_streams, case_info, opt_global)
%% This function will perform clustering with a user-selected algorithm and
% number of clusters. The mass flowrates and the volumes of the reactors are computed 
%
%   *INPUT*:
%       _k_ = number of reactors
%
%       inlet_streams = cell array vector containing the information of the
%                       inlet conditions. The length is equal to the number
%                       of inlets, each cell contain a stream object which
%                       is a struct type array with the following
%                       attributes:
%
%                       stream.id = ID of the stream (BC) from Fluent
%                       stream.T  = Temperature [K]
%                       stream.Y  = mass fractions in the form {'H2:0.5'..}
%                       stream.P  = Pressure [Pa]
%                       stream.Mf = Mass flowrate [kg/s]
%                                            
%       case_info: struct array containing the following fields:
%                   
%                   case_info.Tout = Estimated outlet temperature in K
%                   case_info.Qdot = Estimated heat losses of the system
%
%   *OUTPUT*:
%       _mass flowrates_ = k x k matrix with mass flowrates
%                        across clusters. mass_flowrates(i,j) = mass from i
%                        to j
%
%       _net volumes_ = k x 1 vector of indiivdual reactor volumes
%
%       _bc_mass_flowrarates = k x 2 vector. Element (i,1) is the external
%       mass input in reactor i, element (i,2) is the output from reactor i
%
%       _idx_ = n_cells x 1 vector with the clustering index
%
%       _X_data_ = n_cells x n_vars data matrix (not scaled, not centered)
%                       
%
%       _mass imbalance_ = mass imbalance of the reactors
%
%
%
%   *INSTRUCTIONS*:
%
% In the working folder, you should have a list of .csv files exported
% directly from Fluent. Those files are:
%
%   data_solution: exported variables according to the user, under the name
%   of "data_solution"
%
%   data_volumes: exported cell volumes from Fluent under the name of "data_volumes"
%
%   data_density: exported cell density from Fluent under the name of
%                 "data_density"
%
%   neighbor_cell_flowrates: exported connectivity info and mass flowrates
%                            across neighbor cell using the udf export_neighbor_faces.c
%                            under the name of "neighbor_cell_flowrates"
%
%   Boundary_cells: exported cells at yhe boundary conditions, with
%                   identified BC index and the associated mass flowrate
%
% *AVAILABLE CLUSTERING ALGORITHMS*:
%       _k-means_
%       _local pca_
%       _gaussian mixture models_
%
% *BEFORE USING THE CODE*
%
% You should set-up some important variables for clustering preferences
%
%   select_data: string format, 'all' for all variables, 'major' for major
%   species
%
% The code will detect the zones in the combustor. The nodes which do not
% communicate with each other will be reassigned based on graph algorithms
%
% Each reactor will have several attributes, in the form of a struct array:
%
% R.id = Reactor identifier (integer from 0 to number of reactors)
% R.Type = 'cstr' or 'pfr' (string)
% R.V = volume [cm3]
% R.T = Temperature (1500 K by default if no other information is given)
% R.D = Diameter [m] pfr only
% R.L = Length [m]   pfr only
% R.A = Heat transfer surface [m2]
% R.Tau = Residence time [s]
% R.Qdot = Heat losses [W]

init_time = datetime;

%% Data import section

% Check existence of option DataPath
if isfield(opt_global, 'DataPath') == false
    error('Path to data not available in the options. Impossible to retrieve data');
end

% Import other data
data_volumes  = importdata(append(opt_global.DataPath, 'data_volumes'));
data_density  = importdata(append(opt_global.DataPath, 'data_density'))  ;
data_mass     = importdata(append(opt_global.DataPath, 'Neighbours_cell_flows'));
data_velocity = importdata(append(opt_global.DataPath, 'data_velocity'));
data_solution = importdata(append(opt_global.DataPath, 'data_solution'));

if isfile(append(opt_global.DataPath, 'data_tau')) == true
    disp('Importing also residence time data');
    data_tau = importdata(append(opt_global.DataPath, 'data_tau'));
end

% Check for diffusion data
if isfile(append(opt_global.DataPath, 'data_viscosity'))
    data_viscosity = importdata(append(opt_global.DataPath, 'data_viscosity'));
end

% Extract the data
if isfield(opt_global, 'OptData') == false
    opt_global.OptData = input('Select set of data to import, available options are: all, major, reduced_set: ', 's');
end

% Data matrix (refer to import_data_cfd.m)
[X, labels, dim, coord] = import_data_cfd(opt_global.OptData);
fprintf('Data matrix X is of size %d rows and %d columns \n', size(X, 1), size(X, 2));
start = dim+2; % Columns to start importing the data
            
% Volume data
data_volumes = data_volumes.data;
vol_data = data_volumes(:, start);

% Density data
data_density = data_density.data;
density_data = data_density(:,start);

% Velocity data
vel_data = data_velocity.data;
vel_data = vel_data(:,start);

% Viscosity data
if isfile('data_viscosity')
    viscosity = data_viscosity.data(:,start);
end

% Get mixture fraction
if isfield(opt_global, 'Precond') == true
    load(append(opt_global.DataPath, 'case_info.mat'));
    data_comp = importdata('data_solution');
    comp = data_comp.data;
    comp = comp(:,start+1:end);
    sp_labels = labels(start+1:end);
 
    sp_labels_c = rewrite_fluent_labels(sp_labels);
    Y = mole_to_mass(comp, sp_labels_c);
    f = mixture_fraction(Y, sp_labels, fuel);
    opt_global.f = f;
end

%% Choose algorithm
if isfield(opt_global, 'Alg') == false
    opt_global.Alg = input('Enter clustering algorithm: k-means or lpca: ', 's');
end
alg = opt_global.Alg;

% Center the data
if isfield(opt_global, 'Center') == false
    cent = input('Enter 1 to center the data, 0 to proceed without centering: ');
else
    cent = opt_global.Center;
end

% Select scaling criteria
if isfield(opt_global, 'Scale') == false
    scal = input('Scaling criteria, available choises are: no, auto, range, pareto, vast, level, max, : ', 's');
else
    scal = opt_global.Scale;
end
if strcmp(scal, 'auto') == true
    scal_crit = 1;
elseif strcmp(scal, 'range') == true
    scal_crit = 2;
elseif strcmp(scal, 'pareto') == true
    scal_crit = 3;
elseif strcmp(scal, 'vast') == true
    scal_crit = 4;
elseif strcmp(scal, 'level') == true
    scal_crit = 5;
elseif strcmp(scal, 'max') == true
    scal_crit = 6;
elseif strcmp(scal, 'no') == true
    scal_crit = 0;
else
    disp('No criteria has been selected, auto scale by default');
    scal_crit = 1;
end

% Center and scale
X_center = center(X, cent);
X_scaled = scale(X_center, X, scal_crit);

time = date;
str = append(alg, '_', num2str(k), 'clust_', scal, '_', time);

% Select options for clustering
if strcmp(alg, 'lpca') == true

    if isfield(opt_global, 'StopRule') == false
        opt_global.StopRule = input('Enter stopping rule for lpca. Available are var, eigs: ', 's');
    end

    if strcmp(opt_global.StopRule, 'var') == true
        if isfield(opt_global, 'Inputs') == false
            opt_global.Inputs = input('Select the amount of variance to retain: ');
        end

    elseif strcmp(opt_global.StopRule, 'eigs') == true
        if isfield(opt_global, 'Inputs') == false
            opt_global.Inputs = input('Select number of eigenvectors to retain: ');
        end
    end
end

% Perform clustering
disp('Performing unsupervised clustering...');
clust_time_start = datetime;

idx = custom_cluster(k, X, labels, alg, opt_global);
scatter(coord(:,2), coord(:,1), 5, idx, 'filled');
pause;
close all;

clust_time_end   = datetime;
infos.ClusteringTime = clust_time_end - clust_time_start;

disp('Clustering completed, computing volumes...');

%% Connectivity study

% Remove the boundary cells from the connectivity matrix (refers to remove_boundary.m)
data_mass = data_mass.data;
data_mass = data_mass(:,4:end);
mass_data = remove_boundary(data_mass);

% Build the graph (refers to create_graph.m)
nodes = [1:1:length(X)]';
conn_matrix = mass_data(:,1:2);
G = create_graph(nodes, conn_matrix);

% Analyze the connectivity (refers to check_connectivity.m and clustering.m)
% Partition the cell id's
cell_clust = clustering(nodes, idx);

% Check for empty clusters
idx_new = idx;
for j = 1 : length(cell_clust)
    if isempty(cell_clust{j})
        mess = append('Cluster ', num2str(j), ' is empty');
        disp(mess);
        for i = 1 : length(idx)
            if idx_new(i) > j
                idx_new(i) = idx_new(i) - 1;
            end
        end
    end
end

% Evaluate the new number of clusters
k_new = max(idx_new);
idx = idx_new;
mess = append('Now we have ', num2str(k_new), ' clusters');
disp(mess);

% Re-group the cells
cell_clust = clustering(nodes, idx_new);

% Check connectivity
[all_connected, n_sub_graph] = check_connectivity(G, cell_clust);

disp('Connectivity study completed. Re-assigning nodes...');

% Re-assign nodes for connectivity
graph_time_start = datetime;

% Re-assign options
opt_r.Criterion = 'Euclidean';
opt_r.X = X_scaled;
opt_r.ReassignCriterion = 'volume';
opt_r.CellsVolume = vol_data;
opt_r.VolumeThreshold = 0.01;

idx_new = reassign_nodes(nodes, G, idx, opt_r);
idx = idx_new;
k_new = max(idx_new);
if k_new ~= k
    clust_removed = k - k_new;
    mess = append(num2str(clust_removed), ' has been removed');
    k = k_new;
end
% Get infos
graph_time_end = datetime;
infos.GraphTime = graph_time_end - graph_time_start;

%% Calculate volumes and mass of clusters
% Initialize volumes and mass of reactors
net_volumes = zeros(k,1);
mass = zeros(k,1);

% If 2d-axisimmetric volumes must be corrected
if size(coord,2) == 2
     vol_data = vol_data; % * 2 * pi;
end

% Partititon the data and compute volumes and mass
clust_vol = clustering(vol_data, idx);
clust_dens = clustering(density_data, idx);

% Sum the volume of the cell in each cluster
for j = 1 : k
    net_volumes(j) = sum(clust_vol{j})*1e6;                 % cm3
    mass(j) = 1e-6*net_volumes(j)*mean(clust_dens{j});      % kg
end

%% Remove the clusters below 0.01% of the total volume
opt_global.RemoveSmall = true;
opt_global.RemoveThreshold = 1e-6;
if isfield(opt_global, 'RemoveSmall') == true
    fprintf('Removing small clusters... \n');

    if isfield(opt_global, 'RemoveThreshold') == true
        f_vol = opt_global.RemoveThreshold;
    else
        f_vol = input('Enter the fraction threshold for removing small clusters (default 1e-6): ');
    end

    vol_thresh = sum(net_volumes) * f_vol;
    fprintf('\n Removing clusters smaller than %d cm3... \n', vol_thresh);
    [idx_new, new_volumes, n_small] = remove_small_clusters(net_volumes, vol_thresh, idx, G);
    infos.NSmall = n_small;
    
    k_new = max(idx_new);
    k = k_new;
    
    % Update variables
    idx = idx_new;
    net_volumes = new_volumes;
    
    fprintf('\n Nodes reassigned \n');
    fprintf('\n Number of clusters obtained: %d \n \n', k_new);
    
    % Check the connectivity again
    cell_clust = clustering(nodes, idx);
    [all_connected, n_sub_graph] = check_connectivity(G, cell_clust);
end

%% Compute mass flowrates across reactors
disp('Calculating mass flowrates...');

% Initialize mass flowrate matrix. Mass flowrate (i,j) is the mass flowing
% from reactor i to j
mass_flowrates = zeros(k,k);

% Scanning through neighbor cells
for j = 1 : size(mass_data, 1)
    id_1 = mass_data(j,1);
    id_2 = mass_data(j,2);
    
    % If neighbor cells belong to different cluster, compute the mass
    % flowrate
    if id_2 ~= 0 && id_2 ~= -1 && isnan(id_1) == false
        clust_1 = idx(id_1);
        clust_2 = idx(id_2);
        if clust_1 ~= clust_2
            mi = mass_data(j,3);
            
            % If mass flowrate > 0, the flow is going from clust_1 to
            % clust_2, so store it in mass_flowrates(clust_1, clust_2); if
            % the mass flowrate < 0 is the opposite, store it in
            % mass_flowrates(clust_2, clust_1);
            
            if mi > 0
                mass_flowrates(clust_1, clust_2) = mass_flowrates(clust_1, clust_2) + mi;
            else
                mass_flowrates(clust_2, clust_1) = mass_flowrates(clust_2, clust_1) - mi;
            end
        end
    end
end

%% Display results
% Plot the results
if isfield(opt_global, 'Plot')
    plt = opt_global.Plot;
else
    plt = input('Plot the clustering results? Enter true or false: ');
end

if plt == true
    output = plot_clustering(coord, idx, dim);
end

%% Correcting mass flowrates for mass balance check

% Calculate the mass flowrate at the boundaries (inlets and outlets)
data_boundaries = importdata('Boundary_cells');

% Calculate the inlets and mass flowrates at the boundaries
[inlets, bc_mass_flowrates] = calc_inlets(idx, inlet_streams, data_boundaries);

% Calculate reactors areas
if opt_global.SolveEnergy == true
    Ar = get_face_area(data_boundaries, idx);
end
        
% Evaluate both global and local mass unbalance
disp('Checking mass balance...');

% Select true to correct the mass flowrates, false to keep the original
% calculated value
if isfield(opt_global, 'MassAdjust') == false
    correction = false;
    fprintf('Options for adjusting mass flowrates not specified. They will not be modified \n');
else 
    correction = opt_global.MassAdjust;
end

switch correction
    case true
        alpha = mass_split(mass_flowrates, bc_mass_flowrates);
        for i = 1 : k
            for j = 1 : k
                if alpha(i,j) < 1e-3
                    alpha(i,j) = 0;
                end
            end
        end
        

        % Solve the system to get the corrected mass flowrates
        A  = eye(k) - alpha';
        b  = bc_mass_flowrates(:,1);
        mf = linsolve(A,b);
        mass_flowrates = alpha.*mf;     
end

% Check if there are clusters with exchanging no mass
k_new = k;
ind = [1:1:k]';

for i = 1 : k
        if isempty(find(mass_flowrates(i,:) ~= 0)) && isempty(find(mass_flowrates(:,i) ~= 0))
            % Remove row and column
            removerows(ind, 'ind', i);
            k_new = k - 1;
            fprintf('No mass exchange has been detected in cluster %d, it will be deleted', i);
        end
end

mass_flowrates = mass_flowrates(ind, ind);
k = k_new;
            
[global_imbalance, local_imbalance] = net_balance(mass_flowrates, bc_mass_flowrates);
fprintf('Global relative imbalance in the network = %E \n', global_imbalance);
fprintf('Relative mass imbalance in each cluster: \n');
for j = 1 : k
    fprintf('%E \t', local_imbalance(j));
end
fprintf('\n');

mass_corrected = mass_flowrates;

%% Save the results into a newly created folder
dir_title = str;
mkdir(dir_title);
cd(dir_title);
infos.MainDir = dir_title;

save clustering_results
if plt == true
    figure(1); hold on; saveas(gcf, 'Cluster_output.png');
end

% Save the results
save clustering_results

% Write output file as .txt
if isfield(opt_global, 'WriteFile') == true
    save idx.txt idx -ascii;
    save X.txt X -ascii;
end

%% Write out the files for NetSMOKE++
if isfield(opt_global, 'KineticMech') == true
    netsmoke_dir = append(dir_title, '_', opt_global.KineticMech);
else
    netsmoke_dir = append(dir_title, '_gri3.0');
end

infos.NetSmokeDir = netsmoke_dir;

% Create a new folder
mkdir(netsmoke_dir);
cd(netsmoke_dir);

val = data_solution.data;
T = val(:,start);

% Data necessary for NetSmoke
% Calculate the total outflow
bc_out = find(bc_mass_flowrates(:,2) < 0);
m_out = -sum(bc_mass_flowrates(:,2));

% Calculate the inlets and mass flowrates at the boundaries
[inlets, bc_mass_flowrates] = calc_inlets(idx, inlet_streams, data_boundaries);

% Group coordinates and temperature of the clusters
coord_clust = clustering(coord, idx);
T_clust = clustering(T, idx);

% Calculate the mass of each cell and multiply it by the temperature in
% order to get the mass-weighted average at the end
mass_cell = vol_data.*density_data;
weight_T = T.*mass_cell;
weight_clust = clustering(mass_cell, idx);
weight_T_clust = clustering(weight_T, idx);

% For each reactor create an object with the required attributes and a
% global list of all the reactors
R_list = cell(k,1);

% Energy equation to solve?
if isfield(opt_global, 'SolveEnergy') == false
    en_eq = input('Solve the energy equation? enter true or false: ');
else
    en_eq = opt_global.SolveEnergy;
end

a = 1000;   % Multiplicative factor for mass flowrates

x_clust = clustering(coord(:,1),idx);

% Get variance of temperature in the clusters
if opt_global.DataVariance == true
    data_variance = importdata('../../data_Trms');
    Tvar = data_variance.data(:,start);
    Tvar_clust = clustering(Tvar, idx);
    fprintf('Temperature variance data are present and they are being processed... \n');
end

%% Check for diffusive fluxes between reactors
write_input_diffusion = false;
if isfield(opt_global, 'Diffusion')
    % Calculate distance
    DD = calc_cluster_distance(coord, idx);
    % Subgraphs
    [~, ~, ~, ~, ~, subgraphs] = ...
        connectivity_study(G, val(:,1), idx, false, coord);
    % Calculate contact surface
    data_neighbors = importdata('../../Neighbours_cell_flows');
    val_neighbors  = data_neighbors.data;
    Ab = calc_boundary_surface(G, subgraphs, idx, val(:,1), val_neighbors);
    Ab = Ab * 2 * pi;
    % Calculate diffusion fluxes
    Y = data_solution.data(:, start+1:end);
    Dm = estimate_diffusion_fluxes(Y, DD, Ab, viscosity, idx);
    disp('Diffusive fluxes successfully estimated');
    % Flag to write input for diffusion
    write_input_diffusion = true;
end

%% Create the reactors objects
switch en_eq
    case true

        for i = 1 : k
            R.Type = 'cstr';            % Reactor type
            R.id = i-1;                 % Reactor number identifier (from 0)
            
            % If it is a CSTR then we need volume and mass flowrate
            if strcmp(R.Type, 'cstr') == true
                R.V = net_volumes(i);
                R.Mf = sum(mass_flowrates(i,:)) + bc_mass_flowrates(i);
                R.P = case_info.P;
                R.T = sum(weight_T_clust{i})/sum(weight_clust{i});

                if max(x_clust{i}) < 0.0
                    R.isothermal = true;
                else
                    if R.T > 1000.0 && opt_global.KineticCorrections == true
                        R.KineticCorrections = true;
                        R.Tmean = R.T;
                        % Check if data of temperature variance is present
                        if opt_global.DataVariance == false
                            R.Tvar  = opt_global.RelativeFluctuations * R.Tmean;
                        else
                            R.Tvar = mean(Tvar_clust{i});
                        end
                    end
                end

                % Check if option to initialize the solution from CFD is
                % active
                if isfield(opt_global, 'InitializeComposition') == true
                    % Extract species
                    val = data_solution.data;
                    % Select major species
                    mj = {'CH4', 'O2', 'H2O', 'H2', 'OH'};
                    id = zeros(length(mj),1);
                    for s = 1 : length(mj)
                        id(s) = find(strcmp(sp_labels_c, mj(s)));
                    end
                    cc = comp(:,id);
                    cc_clust = clustering(cc, idx);
                    cc_mean = mean(cc_clust{i});
                    Ys = cell(length(mj),1);
                    for s = 1 : length(cc_mean)
                        Ys{s} = append(mj{s}, ':', num2str(cc_mean(s)));
                    end
                    R.Y = Ys;
                    R.basis = 'mole';
                else   
                    R.Y = {'O2:0.21', 'N2:0.79'};   % Mole fractions
                    R.basis = 'mole';
                end
                
                % If it's an inlet, keep it isothermal
                if isempty(inlets{i}) == false
                    R.isothermal = true;
                    R.V = 0.01;

                % Check if the reactor is very small, in that case change
                % the ODE settings
                else
                    if R.V < 10 %cm3
                        R.odesettings = true;
                        R.reltol = 1e-12;
                        R.abstol = 1e-12;
                    end
                end

                % Check if lateral area of reactor is more than 1% of total
                % area, in that case reactor(i) has thermal exchange
                if Ar(i)/sum(Ar) > 0.01
                    th_exch = false;
                    switch th_exch
                        case true
                            R.Qdot = case_info.Qdot*Ar(i)/sum(Ar);
                            R.Tout = case_info.Tout;
                            R.Tenv = 300;
                            R.A = 2.44*Ar(i)/sum(Ar);
                        case false
                            R.isothermal = true;
                    end
                end
                
            % If it is a PFR we need the mass flowrate, residence time and diameter
            else
                R.Mf  = sum(mass_flowrates(i,:)) + bc_mass_flowrates(i);                % kg/s
                R.Tau = mass(i)/R.Mf;                                                   % s
                R.L = abs(max(coord_clust{i}(:,1)) - min(coord_clust{i}(:,1)))*1000;    % mm
                R.P = case_info.P;
                R.T = 1500;
                R.Y = {'O2:0.21', 'N2:0.79'};
            end
            
            output = write_reactor_input(R);
            
            % Append R to the list
            R_list{i} = R;
            
            clear R;
        end

    case false

        for i = 1 : k
            R.Type = 'cstr';            % Reactor type
            R.id = i-1;                 % Reactor number identifier (from 0)
            
            % If it is a CSTR then we need volume and mass flowrate
            if strcmp(R.Type, 'cstr') == true

                % Set main attributes
                R.V = net_volumes(i);
                R.Mf = sum(mass_flowrates(i,:)) + bc_mass_flowrates(i);
                R.P = case_info.P;

%                 R.T = sum(weight_T_clust{i})/sum(weight_clust{i});
                R.T = max(T_clust{i});

                R.Y = {'O2:0.21', 'N2:0.79'};   % Mole fractions                
                R.isothermal = true;

                % If it's an inlet, keep it isothermal
                if isempty(inlets{i}) == false
                    R.V = 0.01;

                % Check if the reactor is very small, in that case change
                % the ODE settings
                else
                    if R.V < 100 %cm3
                        R.odesettings = true;
                        R.reltol = 1e-12;
                        R.abstol = 1e-12;
                    end
                end

                % Check for kinetic corrections
                if R.T > 1000 && opt_global.KineticCorrections == true
                    R.KineticCorrections = true;
                    R.Tmean = R.T;
                    % Check if data of temperature variance is present
                    if opt_global.DataVariance == false
                        R.Tvar  = opt_global.RelativeFluctuations * R.Tmean;
                    else
                        R.Tvar = max(Tvar_clust{i});
                    end
                end

                % Write input of the single reactor
                if write_input_diffusion == true
                    R.TurbDiff = mean(viscosity(idx==i));
                end
                
            % If it is a PFR we need the mass flowrate, residence time and diameter
            else
                R.Mf  = sum(mass_flowrates(i,:)) + bc_mass_flowrates(i);            % kg/s
                R.Tau = mass(i)/R.Mf;                                                     % s
                R.L = abs(max(coord_clust{i}(:,1)) - min(coord_clust{i}(:,1)))*1000;     % mm
                R.P = case_info.P;
                R.T = mean(T_clust{i});
                R.Y = {'O2:0.21', 'N2:0.79'};
            end
            
            % Write input of the single reactor
            output = write_reactor_input(R);
            
            % Append R to the list
            R_list{i} = R;
            
            clear R;
        end

end

% Create the inlets reactors
n_in = length(find(bc_mass_flowrates(:,1) ~= 0));   % Number of inlets

new_connections = zeros(n_in, 3);
count = 0;
for i = 1 : k
    if bc_mass_flowrates(i,1) ~= 0
        % Update counter
        count = count + 1;
        
        % Then create a "fake" reactor
        R.Type = 'cstr';
        R.id = k + count - 1;
        R.V = 1;                            % cm3
        R.Mf = bc_mass_flowrates(i,1);      % kg/s
        R.Y  = inlets{i}.Y;                 % mass fractions
        R.T  = inlets{i}.T;                 % K
        R.P  = inlets{i}.P;                 % Pa
        R_list{k+count} = R;
        
        % Write input
        output = write_reactor_input(R);
        
        % Update new connections
        new_connections(count,1) = k + count;
        new_connections(count,2) = i;
        new_connections(count,3) = R.Mf;

    end
end

% Create the new mass flowrate matrix
mass_netsmoke = zeros(k+n_in);
for i = 1 : k
    for j = 1 : k
        mass_netsmoke(i,j) = mass_flowrates(i,j)*a;
    end
end

% Then add the mass flowrates from the fake reactors
for i = 1 : n_in
    mass_netsmoke(new_connections(i,1), new_connections(i,2)) = new_connections(i,3)*a;
end

% Create the new bc_mass_flowrates vector
m_ext = zeros(k+n_in, 2);
for i = 1 : k
    m_ext(i,2) = -bc_mass_flowrates(i,2)*a;
end

for i = 1 : n_in
    m_ext(new_connections(i,1),1) = new_connections(i,3)*a;
end

% If diffusion is active we need to adjust dimensions of areas and
% distances
if isfield(opt_global, 'Diffusion')
    Anew = zeros(size(mass_netsmoke));
    Dnew = ones(size(mass_netsmoke));

    [sz1, sz2] = size(Ab);
    Anew(1:sz1, 1:sz2) = Ab;
    Dnew(1:sz1, 1:sz2) = DD;
end

%% Write NetSMOKE++ inputs
if isfield(opt_global, 'KineticMech') == true
    kin_mech = opt_global.KineticMech;
else
    kin_mech = 'gri3.0';
end

switch kin_mech
    case 'SanDiego'
        opt_out.Chemfile = '/Users/matteosavarese/Desktop/Dottorato/Kinetics/San_Diego_NOx/SanDiego.CKI';
        opt_out.Thermofile = '/Users/matteosavarese/Desktop/Dottorato/Kinetics/San_Diego_NOx/thermo.dat';
        opt_out.KineticPath = true;
    case 'Polimi'
        opt_out.Chemfile = '/Users/matteosavarese/Desktop/Dottorato/Kinetics/Polimi_mech/Polimi_c1c3_nox_red.CKI';
        opt_out.Thermofile = '/Users/matteosavarese/Desktop/Dottorato/Kinetics/Polimi_mech/thermo_c1c3_nox.dat';
        opt_out.KineticPath = true;
    case 'gri3.0'
        opt_out.Chemfile = '/Users/matteosavarese/Desktop/Dottorato/Kinetics/GRI3.0/gri30.CKI';
        opt_out.Thermofile = '/Users/matteosavarese/Desktop/Dottorato/Kinetics/GRI3.0/thermo.dat';
        opt_out.KineticPath = true;
    case 'gri2.11'
        opt_out.Chemfile = '/Users/matteosavarese/Desktop/Dottorato/Kinetics/GRI_2_11/Kinetics.cki.txt';
        opt_out.Thermofile = '/Users/matteosavarese/Desktop/Dottorato/Kinetics/GRI_2_11/Thermo.dat.txt';
        opt_out.KineticPath = true;
    case 'Stagni'
        opt_out.Chemfile = '/Users/matteosavarese/Desktop/Dottorato/Ammonia_SNCR_Project/Open_Smoke_simulations/Stagni/chem.dat';
        opt_out.Thermofile = '/Users/matteosavarese/Desktop/Dottorato/Ammonia_SNCR_Project/Open_Smoke_simulations/Stagni/therm.dat';
        opt_out.KineticPath = true;
    case 'Otomo'
        opt_out.Chemfile = '/Users/matteosavarese/Desktop/Dottorato/Ammonia_SNCR_Project/Open_Smoke_simulations/Otomo/otomo.inp';
        opt_out.Thermofile = '/Users/matteosavarese/Desktop/Dottorato/Ammonia_SNCR_Project/Open_Smoke_simulations/Otomo/otomo_therm.dat';
        opt_out.KineticPath = true;
end

opt_out.SpeciesMonitor = {'NO', 'NH3'};

% Check if diffusion input must be written
if write_input_diffusion == false
    output = write_global_input(R_list, mass_netsmoke, m_ext, opt_out);
else
    output = write_global_input_diffusion(R_list, mass_netsmoke, m_ext, Anew, Dnew, opt_out);
end

% Write .sh file
if isfile('Run.sh') == false
    bat = fopen('Run.sh', 'w');
    fprintf(bat, '/Users/matteosavarese/Desktop/Dottorato/OpenSmoke/Compilazione_KPP/NetSmoke_KinCorr-Master/SeReNetSMOKEpp-master/projects/Linux/SeReNetSMOKEpp.sh --input input.dic');
    fclose(bat);
end

% Run the simulation
if isfield(opt_global, 'RunSimulation') == false
    run_sim = input('Run the simulation? Enter true or false: ');
else
    run_sim = opt_global.RunSimulation;
end

if run_sim == true
    net_time_start = datetime;
    ! sh Run.sh
    net_time_end   = datetime;
    infos.SimTime = net_time_end - net_time_start;
end

disp('Global input written');

%% Measure time for running
end_time = datetime;
dur_time = end_time - init_time;
infos.TotalTime = dur_time;

fprintf('Unsupervised clustering time: ');
disp(infos.ClusteringTime);

fprintf('Graph reassignment time: ');
disp(infos.GraphTime);

if run_sim == true
    fprintf('Simulation time: ');
    disp(infos.SimTime);
end

fprintf('Total execution time: ');
disp(infos.TotalTime);

%% Post-processing
if isfield(opt_global, 'PostProcessing') == true
    pp = opt_global.PostProcessing;
else
    pp = input('Post processing the results? Enter true or false: ');
end

pp = true;
if pp == true
    opt_global.ExportSpecies = {'NO', 'NO2', 'N2O', 'NH3'};
    output_species = post_process_results('Output', 'input.dic', opt_global);
    fprintf('\n Output species: \n');
    for i = 1 : length(output_species.Species)
        fprintf('%s = %E \n', output_species.SpeciesName{i}, output_species.Species(i));
    end
end

if opt_global.PlotContour == true && opt_global.RunSimulation == true
    [Y_net, Y_net_contour] = compare_contour(data_solution.data(:,4), 'T[K]', 'Output', idx, coord);
    [infos.Y_net_NO] = plot_network_contour(idx, coord, 'NO', 'Output/');
    
    % Plot NO contour
    figure;
    scatter(coord(:,2), coord(:,1), 20, infos.Y_net_NO, 'filled');
    hold on;
    scatter(-coord(:,2), coord(:,1), 20, infos.Y_net_NO, 'filled');
    colormap('hot(20)');
    cb = colorbar;
    cb.Label.String = 'NO [mol.f]';
    fig = gcf; fig.Units = 'centimeters';
    fig.Position = [15 15 18 14];
end

if isfield(opt_global, 'PlotMixtureFraction') == true
    plt_mf = true;
else
    plt_mf = input('Plot mixture fraction vs T? Enter true or false: ');
end

if plt_mf == true

    figure;
    scatter(opt_global.f, T, 10, idx, 'filled');
    cb = colorbar;
    cb.Label.String = 'Cluster index';
    cmap = append('parula(', num2str(k), ')');
    colormap(cmap);
    xlabel('F [-]');
    ylabel('T [K]');
end

opt_global.PhysicalIndex = true;
if isfield(opt_global, 'PhysicalIndex') == true

    % Plot standard deviation, mean and range
    Tmean = zeros(k,1);
    Tstd  = zeros(k,1);
    Trange = zeros(k,1);
    phc = zeros(k,1);
    phc_dev = zeros(k,1);
    for i = 1 : k
        Tmean(i)  = mean(T_clust{i});
        Tstd(i)   = std(T_clust{i});
        Trange(i) = max(T_clust{i}) - min(T_clust{i});
        phc(i) = Trange(i)/Tmean(i);
    end

    figure;
    subplot(2,1,1);
    errorbar([1:1:k]', Tmean, Tstd, 'k', 'LineStyle','none', 'Marker','o', 'MarkerFaceColor','k');
    ax = gca; ax.TickLabelInterpreter = 'latex';
    xlabel('Cluster index');
    ylabel('T [K]');

    % Calculate physical homogeneity coefficient
    subplot(2,1,2);
    plot([1:1:k]', phc, 'LineStyle','none', 'LineWidth',2, 'Marker','o', 'MarkerFaceColor','k', 'Color','k');
    xlabel('Cluster index');
    ylabel('PHC');
    ax = gca; ax.TickLabelInterpreter = 'latex';
end

if isfield(opt_global, 'PlotProfile') == true
    if opt_global.PlotProfile == true
        disp('Plotting network profile...');
        [y_axis, x_net, y_cfd, x_cfd] = plot_network_profile(idx, coord, T, 'Output/', plt);
    end
end

% Plot NO contour
plotNO = false;
if plotNO == true
    Y_net = plot_network_contour(idx, coord, 'NO', 'Output/');

    % Restrict domain to y > 0
    id0 = find(coord(:,1) > 0);
    coord_0 = coord(id0,:);
    Y_net0 = Y_net(id0);

    % Meshgrid
    xmesh = linspace(0, max(coord(:,2)), 500);
    ymesh = linspace(0, max(coord(:,1)), 500);
    [XX, YY] = meshgrid(xmesh, ymesh);

    % Interpolate solution
    Ynet_cont = griddata(coord_0(:,2), coord_0(:,1), Y_net0, XX, YY);

    % Plot contour
    figure;
    contourf(XX, YY, Ynet_cont, 50, 'LineStyle','none');
    hold on;
    contourf(-XX, YY, Ynet_cont, 50, 'LineStyle','none');
    cmap = brewermap(25, '-RdBu');
    colormap(cmap);
    fig = gcf; fig.Units = 'centimeters';
    fig.Position = [15 15 10 18];
    cb = colorbar;
end

% Compare contours
compare_contours = true;
if compare_contours == true
    
    % Select index of the variable to plot
    id = 11;
    varname = 'H';
    Y_cfd = data_solution.data(:,id);
    [Y_net, Y_net_contour] = compare_contour(Y_cfd, varname, 'Output/', idx, coord);

    close;

    % Plot contours
    % Restrict domain to y > 0
    id0 = find(coord(:,1) > 0);
    coord_0 = coord(id0,:);
    Y_net0 = Y_net_contour(id0);
    Y_cfd0 = Y_cfd(id0);

    % Meshgrid
    xmesh = linspace(0, max(coord(:,2)), 500);
    ymesh = linspace(0, max(coord(:,1)), 500);
    [XX, YY] = meshgrid(xmesh, ymesh);

    % Interpolate solution
    Ynet_cont = griddata(coord_0(:,2), coord_0(:,1), Y_net0, XX, YY);
    Ycfd_cont = griddata(-coord_0(:,2), coord_0(:,1), Y_cfd0, -XX, YY);

    % Plot contour
    figure;
    contourf(XX, YY, Ynet_cont, 50, 'LineStyle','none');
    hold on;
    contourf(-XX, YY, Ycfd_cont, 50, 'LineStyle','none');
    cmap = brewermap(50, '-RdBu');
    colormap(cmap);
    fig = gcf; fig.Units = 'centimeters';
    fig.Position = [15 15 14 20];
    cb = colorbar;

    xlabel('x (mm)');
    ylabel('y (mm)');
    cb.Label.String = varname;

    xlim([-0.1 0.1]);
    
    if write_input_diffusion == true
        figname = append('Comparison_', varname, '_diff.png');
        exportgraphics(fig, figname, 'Resolution',600);
    else
        figname = append('Comparison_', varname, '_no_diff.png');
        exportgraphics(fig, figname, 'Resolution',600);
    end


end


save info_simulation infos;

cd ../../

end

