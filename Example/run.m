%% Make sure all the functions are in the path
addpath(genpath('/Users/matteosavarese/Desktop/Dottorato/Github/NetGEN'));

%% Define options
% % LPCA options
% opt_global.Alg = 'k-means';
% opt_global.StopRule = 'eigs';
% opt_global.Inputs = 2;
% % opt_global.Precond = false;
% opt_global.Plot = true;
% opt_global.Initialization = 2; % 1 = random, 2 = uniform, 3 = FPCA

% Simulation options
opt_global.KineticMech = 'gri2.11'; % gri3.0, gri2.11, SanDiego, Polimi, Stagni, Otomo
opt_global.PostProcessing = false;
opt_global.SolveEnergy = false;
opt_global.RunSimulation = true;
opt_global.PlotContour = false;
opt_global.PlotProfile = false;
opt_global.KineticCorrections = true;
% opt_global.RelativeFluctuations = 0.20;
opt_global.DataVariance = true;

% Data pre-processing options
opt_global.OptData = 'reduced_set';
opt_global.Mixing = false;
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

% Diffusion between reactors
% opt_global.Diffusion = true;

% Path to data
opt_global.DataPath = '/Users/matteosavarese/Desktop/Dottorato/Github/NetGEN/Example/';

%% Run the simulation
load case_info.mat
[mass_corrected, net_volumes, bc_mass_flowrates, idx, X, infos] = ...
    network_eval(15, inlet_streams, case_info, opt_global);