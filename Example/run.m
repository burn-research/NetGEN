%% Make sure all the functions are in the path
addpath(genpath('/Users/matteosavarese/Desktop/Dottorato/Github/NetGEN'));

%% Define options

% Data import options
opt_global.VolumeFile = 'data_volumes';
opt_global.DensityFile = 'data_density';
opt_global.ConnectivityFile = 'Neighbours_cell_flows';
opt_global.AngleFile = 'data_angle';
opt_global.TauFile   = 'data_tau';
opt_global.TVarianceFile = 'data_variance';
opt_global.VelocityFile = 'data_velocity';
opt_global.BoundaryCellsFile = 'Boundary_cells';
opt_global.SolutionFile = 'data_solution';
opt_global.ViscosityFile = 'data_viscosity';
opt_global.Basis = 'mol';

% Simulation options
opt_global.KineticMech = 'gri2.11'; % gri3.0, gri2.11, SanDiego, Polimi, Stagni, Otomo
opt_global.PostProcessing = false;
opt_global.SolveEnergy = false;
opt_global.RunSimulation = false;
opt_global.PlotContour = false;
opt_global.PlotProfile = false;
opt_global.KineticCorrections = true;
opt_global.DataVariance = true;

% Data pre-processing options
opt_global.OptData = 'reduced_set';
opt_global.Center = 1;
opt_global.Scale = 'auto';

% K-Means options
opt_global.Alg = 'k-means';
opt_global.Start = 'uniform';
opt_global.MaxIter = 1000;

% Write output file 
opt_global.WriteFile = true;

% Options for reassigning nodes
opt_global.VolumeThreshold = 0.01;
opt_global.ReassignCriterion = 'volume';
% Re-assign options
opt_global.Criterion = 'Euclidean';
opt_global.ReassignCriterion = 'volume';
opt_global.VolumeThreshold = 0.01;

% Diffusion between reactors
opt_global.Diffusion = true;

% Path to data
opt_global.DataPath = '/Users/matteosavarese/Desktop/Dottorato/Github/NetGEN/Example/';

%% Run the simulation
load case_info.mat
[mass_corrected, bc_mass_flowrates, R_list, idx, X, infos] = ...
    network_eval(15, inlet_streams, case_info, opt_global);