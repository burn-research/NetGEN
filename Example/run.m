%% Make sure all the functions are in the path
addpath(genpath('/Users/matteosavarese/Desktop/Dottorato/Github/NetGEN'));

set(0, 'defaultaxesfontsize', 20);
set(0, 'defaulttextfontsize', 20);
set(0, 'defaulttextinterpreter', 'latex');
%% Define options
clear all;

% Path to data
opt_global.DataPath = '/Users/matteosavarese/Desktop/Dottorato/Github/NetGEN/Example/DataCFD/';

% Data import options
opt_global.VolumeFile = 'data_volumes';
opt_global.DensityFile = 'data_density';
opt_global.ConnectivityFile = 'Neighbours_cell_flows';
opt_global.AngleFile = 'data_angle';
opt_global.TauFile   = 'data_tau';
opt_global.TVarianceFile = 'data_Trms';
opt_global.VelocityFile = 'data_velocity';
opt_global.BoundaryCellsFile = 'Boundary_cells';
opt_global.SolutionFile = 'data_solution';
opt_global.Basis = 'mol';

% Simulation options
opt_global.KineticMech = 'SanDiego'; % gri3.0, gri2.11, SanDiego, Polimi, Stagni, Otomo
opt_global.SolveEnergy = true;
opt_global.RunSimulation = false;
opt_global.KineticCorrections = true;
opt_global.DataVariance = true;
opt_global.InitComp = true;

% Data pre-processing options
opt_global.OptData = 'reduced_set';
opt_global.Center = 1;
opt_global.Scale = 'auto';

% K-Means options
opt_global.Alg = 'k-means';
opt_global.Start = 'plus';
opt_global.MaxIter = 1000;

% Write output file 
opt_global.WriteFile = true;

% Options for reassigning nodes
opt_global.VolumeThreshold = 0.001;
opt_global.ReassignCriterion = 'volume';
opt_global.Criterion = 'Euclidean';
% Options for reassigning clusters
opt_global.ClusterReassignment = true;
opt_global.ClusterReassignmentCriterion = 'volumes';

% Reordering reactors
opt_global.ReorderReactors = true;

%% Run the simulation
load case_info.mat

% Select number of clusters
k = 20;     % Number of clusters
ndim = 2;   % If case is 2D or 3D

fname1 = append('k-means_', num2str(k), 'clust_', opt_global.Scale, '_', date);
fname2 = append('k-means_', num2str(k), 'clust_', opt_global.Scale, '_', date, '_', opt_global.KineticMech);
fold2 = append(fname1, '/', fname2, '/Output/');

[mass_corrected, bc_mass_flowrates, R_list, idx, X, infos] = ...
     network_eval(k, ndim, inlet_streams, case_info, opt_global);

%% Plot the clustering results
load(append(fname1, '/', 'clustering_results.mat'))
data_all = import_data_all(opt_global.DataPath, ndim, opt_global);

figure;
scatter(data_all.Coordinates(:,1), data_all.Coordinates(:,2), 10, idx, 'filled');
cb = colorbar;
ax = gca; ax.TickLabelInterpreter = 'latex';
xlabel('x [m]'); ylabel('y [m]');
fig = gcf; fig.Units = 'centimeters';
fig.Position = [15 15 18 12];



