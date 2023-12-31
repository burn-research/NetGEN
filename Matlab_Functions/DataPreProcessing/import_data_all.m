function [data_all] = import_data_all(path_to_data, ndim, opt)
%
% function [data_all] = import_data_all(opt)
%
%
%% Preliminary information
if ndim == 2
    st = 4;
elseif ndim == 3
    st = 5;
else
    error('ndim should be 2 or 3 (geometric dimensions)');
end

%% Basic data required for CRN clustering and modelling

% Those are the default files you need to import

% Import volumes
if isfield(opt, 'VolumeFile')
    data_volumes  = importdata(append(path_to_data, opt.VolumeFile));
else
    error('Volume file name not specified');
end

% Import density
if isfield(opt, 'DensityFile')
    data_density  = importdata(append(path_to_data, opt.DensityFile));
else
    error('Density file name not specified');
end

% Import data mass
if isfield(opt, 'ConnectivityFile')
    data_mass     = importdata(append(path_to_data, opt.ConnectivityFile));
else
    error('Connectivity file name not specified');
end

% Data velocity
if isfield(opt, 'VelocityFile')
    data_velocity = importdata(append(path_to_data, opt.VelocityFile));
else
    error('Velocity file name not specified');
end

% Data solution
if isfield(opt, 'SolutionFile')
    data_solution = importdata(append(path_to_data, opt.SolutionFile));
    labels = data_solution.textdata;
    species_names = labels(st+1:end);
else
    error('Solution file name not specified');
end

% Data Boundaries
if isfield(opt, 'BoundaryCellsFile')
    data_boundaries = importdata(append(path_to_data, opt.BoundaryCellsFile));
else
    error('Boundary cells file not specified');
end

% Create struct array with all the data
data_all.Volumes  = data_volumes.data(:,st);            % Volumes of cells
data_all.Density  = data_density.data(:,st);            % Density of cells
data_all.Velocity = data_velocity.data(:,st);           % Velocity magnitude
data_all.Solution = data_solution.data(:,st:end);       % State space data (T mass fractions)
data_all.Coordinates = data_solution.data(:,2:st-1);    % Geometrical coordinates of centroids

% Those should not be touched
data_all.Connectivity = data_mass;
data_all.Boundaries   = data_boundaries;

% Species labels
data_all.SpeciesLabels = species_names;
data_all.Labels        = labels;

%% Customized option

% Import tau data
if isfield(opt, 'TauFile')
    data_tau =  importdata(append(path_to_data, opt.TauFile));
    data_all.Tau = data_tau.data(:,st);
end

% Import viscosity data
if isfield(opt, 'ViscosityFile')
    data_viscosity = importdata(append(path_to_data, opt.ViscosityFile));
    data_all.Viscosity = data_viscosity.data(:,st);
end

% Angle data
if isfield(opt, 'AngleFile')
    data_angle = importdata(append(path_to_data, opt.AngleFile));
    data_all.Angle = data_angle.data(:,st);
end

% Temperature variance data
if isfield(opt, 'TVarianceFile')
    data_variance = importdata(append(path_to_data, opt.TVarianceFile));
    data_all.Variance = data_variance.data(:,st);
end





end

