function [X, labels, dim, coord] = import_data_cfd(opt)

% This function will import the CFD data from the required files and put
% all the matrices and vectors in the struct type data_output
% INPUTS:
%
%   opt = string (available options are: 'all', 'major', 'coord',
%   'reduced_set')
%
% OUTPUTS:
%   X_data = [mxn] matrix of unscaled CFD data
%   labels = array with labels' variables
%   dim    = int. representing the geometrical dimension of the case
%   coord  = mxdim matrix with the coordinates of the mesh

% Check if 'data_solution' file is in the current folder
if isfile('data_solution') == false
    error('Data solution file not found, please make sure you name your solution file like this');
end

% Import the whole dataset
data    =    importdata('data_solution');
labels  =    data.textdata;
X_data  =    data.data;

% Recognize if case is two or three dimensional
% Recognize if the case is 2D or 3D
if strcmp(labels{4}, '    z-coordinate') == true
    dim = 3;
    start = 5;
    disp('Case is 3-Dimensional');
else
    dim = 2;
    start = 4;
    disp('Case is 2-Dimensional');
end

% Get coordinates
coord = X_data(:,2:start-1);

% Select which species you want to keep for the clustering, all, major or
% hydrogen
select_data = opt;
switch select_data
    case 'all'
        X = X_data(:,start:end);
    case 'major'
        major_sp = {'o2', 'oh', 'h2o', 'ch4', 'h2'};
        id_keep = [];
        count = 1;
        for i = 1 : length(major_sp)
            lab = append('molef-', major_sp{i});
            for j = 1 : length(labels)
                current_label = strrep(labels{j}, ' ', '');
                if strcmp(lab, current_label) == true
                    id_keep(count) = j;
                    count = count + 1;
                end
            end
        end
        
        X = [X_data(:,start) X_data(:,id_keep)];
        
    case 'coord'
        major_sp = {'oh', 'o', 'h', 'h2o'};
        id_keep = [];
        count = 1;
        for i = 1 : length(major_sp)
            lab = append('molef-', major_sp{i});
            for j = 1 : length(labels)
                current_label = strrep(labels{j}, ' ', '');
                if strcmp(lab, current_label) == true
                    id_keep(count) = j;
                    count = count + 1;
                end
            end
        end
        X = [X_data(:, 2:3) X_data(:, id_keep)];

    case 'reduced_set'

        major_sp = {'h2o', 'co2', 'co'};
        id_keep = [];
        count = 1;
        for i = 1 : length(major_sp)
            lab = append('molef-', major_sp{i});
            for j = 1 : length(labels)
                current_label = strrep(labels{j}, ' ', '');
                if strcmp(lab, current_label) == true
                    id_keep(count) = j;
                    count = count + 1;
                end
            end
        end

        % Based only on T, f, tau, Tvar which should be contained in new
        % solution files from fluent
        if isfile('case_info.mat') == true
            load case_info.mat;
        else
            error('case_info.mat not found. Please refer to network_eval.m');
        end

        sp_labels = labels(start+1:end);
        sp_labels_c = rewrite_fluent_labels(sp_labels);   % Species labels
        comp = X_data(:,start+1:end);                   % Mole fractions
        Y = mole_to_mass(comp, sp_labels_c);

        % Define the progress variable
        f = sum(Y,2);
        feq = max(f);
        fnorm = f/feq;

        f = mixture_fraction(Y, sp_labels, fuel);

        if isfield(opt, 'Mixing') == true
            data_mixing = importdata('data_mixing');
            mix_time = data_mixing.data(:,4:end);
        end

        % Check for the existence of residence time distribution data
        if isfile('data_tau') == true && isfile('data_velocity') == true && isfile('data_angle') == true && isfile('data_Trms') == true
            data_tau = importdata('data_tau');
            tau = data_tau.data(:,start);
            tau(tau>1) = 1;

            data_var = importdata('data_Trms');
            Tvar = data_var.data(:,4);

            data_velocity = importdata('data_velocity');
            vel = data_velocity.data(:,start);

            data_angle = importdata('data_angle');
            angle = data_angle.data(:,start);
            X = [-log10(tau) X_data(:,start)/10 f vel angle];

            if isfield(opt, 'Mixing') == true
                X = [X mix_time];
            end

        else
            warning('data_tau not found. Rename the file or insert residence time distribution data');
            X = [X_data(:,start) f coord(:,2)];
        end

    case 'pure_h2'
        major_sp = {'o2', 'oh', 'h2o', 'h2', 'h', 'o', 'h2o2', 'ho2', 'n2'};
        id_keep = [];
        count = 1;
        for i = 1 : length(major_sp)
            lab = append('molef-', major_sp{i});
            for j = 1 : length(labels)
                current_label = strrep(labels{j}, ' ', '');
                if strcmp(lab, current_label) == true
                    id_keep(count) = j;
                    count = count + 1;
                end
            end
        end
        
        X = [X_data(:,start) X_data(:,id_keep)];
end

