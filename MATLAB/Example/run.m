%% Make sure all the functions are in the path
addpath(genpath('../../../NetGEN'));

set(0, 'defaultaxesfontsize', 20);
set(0, 'defaulttextfontsize', 20);
set(0, 'defaulttextinterpreter', 'latex');
%% Define options
clear all;

% Path to data
opt_global.DataPath = '../../FluentData/ULBFurnace/';

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
opt_global.SolveEnergy = false;
opt_global.RunSimulation = true;
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
k = 10;     % Number of clusters
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

%% Plot simulation contours

% Select variable to plot
var_name = "OH";
labels = data_all.Labels(4:end);
Y_cfd = data_all.Solution(:,10);
path_to_output = fold2;

[Y_net, Y_net_contour] = compare_contour(Y_cfd, var_name, path_to_output, idx, coord);
cd ../../

% Exclude injector coordinates
id_no_inj = find(coord(:,1)>0);
coord_no_inj = coord(coord(:,1)>0, :);

Yn = Y_net_contour(id_no_inj);
Yc = Y_cfd(id_no_inj);

% Create the mesh
xmin = min(coord_no_inj(:,1));
xmax = max(coord_no_inj(:,1));
ymin = min(coord_no_inj(:,2));
ymax = max(coord_no_inj(:,2));

nps = 250;
xlin = linspace(xmin, xmax, nps);
ylin = linspace(ymin, ymax, nps);

[xq, yq] = meshgrid(xlin, ylin);

% Interpolate data
Yn_cc = griddata(coord_no_inj(:,1), coord_no_inj(:,2), Yn, xq, yq);
Yc_cc = griddata(coord_no_inj(:,1), coord_no_inj(:,2), Yc, xq, yq);

cmap = brewermap(100, "-RdBu");
figure;
% CRN contour
c1 = contourf(yq, xq, Yn_cc, 100, "LineStyle","none");
colormap(cmap)
hold on;
% CFD contour
c2 = contourf(-yq, xq, Yc_cc, 100, "LineStyle","none");
% Middle line
pl = plot([0 0.0], [0.0 0.7], "w--", "LineWidth",1);
% Axis labels
xlabel('x (m)'); ylabel('y (m)');
ax = gca; ax.TickLabelInterpreter = "latex";
% Text
text(-0.30, 0.6, "(CFD)");
text(0.15, 0.6, "(CRN)");
% Colorbar
cb = colorbar;
cb.Label.String = var_name;
cb.Label.Interpreter = "latex";
cb.TickLabelInterpreter = "latex";
% Adjust dimensions
fig = gcf;
fig.Units = 'centimeters';
fig.Position = [12 12 18 20];
% Save
exportgraphics(gca, 'Contour.png', 'Resolution', 600);







