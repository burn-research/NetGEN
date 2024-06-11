function [X, labels, dim, coord] = create_clustering_matrix(data_all, opt)

% This function will import the CFD data from the required files and put
% all the matrices and vectors in the struct type data_output
% INPUTS:
%
%   opt = 's' available options are: 'all', 'major', 'coord', 'reduced_set'
%
% OUTPUTS:
%   X_data = [mxn] matrix of unscaled CFD data
%   labels = array with labels' variables
%   dim    = int. representing the geometrical dimension of the case
%   coord  = mxdim matrix with the coordinates of the mesh

% Import the whole dataset
data    =    data_all.Solution;
labels  =    data_all.Labels;
X_data  =    data;

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
select_data = opt.OptData;

% Check for options
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
        major_sp = {'oh', 'o', 'h'};
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

        major_sp = {'oh', 'h', 'o', 'co2', 'co', 'no'};
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

        sp_labels = data_all.SpeciesLabels;
        sp_labels_c = rewrite_fluent_labels(sp_labels);     % Rewrite species labels
        comp = X_data(:,2:end);                             % Mole or mass fractions

        if isfield(opt, 'Basis') == false
            opt.Basis = input('Mol or mass basis not specified, please specify: ', 's');
        end

        % Transform to mass if necessary
        if strcmp(opt.Basis, 'mol') || strcmp(opt.Basis, 'mole')
            Y = mole_to_mass(comp, sp_labels_c);
        end

        % Define the progress variable
        f = Y(:,1) + Y(:,2);            % Here PV is defined as H2O + CO2
        feq = max(f);
        fnorm = f/feq;

        X = [data_all.Solution(:,1) Y(:,id_keep) fnorm];

        % Check for other fields in data_all
        % Check for tau
        if isfield(data_all, 'Tau')
            tau = data_all.Tau;
             tau = log(tau);       % Apply nonlinear transformation
              tau(tau>1) = 1;

            % Update X
            X = [X tau.^2];
            fprintf('Tau has been added to X \n');
        end

        % Check for velocity
        if isfield(data_all, 'Velocity')
            vel = data_all.Velocity;
            X = [X vel.^1.5];
            fprintf('Velocity has been added to X \n');
        end

        % Check for temperature variance
        if isfield(data_all, 'Variance')
            Tvar = data_all.Variance;
            X = [X Tvar];
            fprintf('Temperature variance has been added to X \n');
        end

%         % Check for velocity angle
%         if isfield(data_all, 'Angle')
%             angle = data_all.Angle;
%             X = [X angle];
%             fprintf('Angle has been added to X \n');
%         end   

    case 'velocity_only'

        vel = data_all.Velocity;
        angle = data_all.Angle;

        X = [vel angle];

        labels = {'velocity', 'angle'};

        if isfield(data_all, 'Tau')
            tau = data_all.Tau;
            tau(tau>1) = 1;         % Clip to max
            tau = -log10(tau);      % Apply nonlinear transformation

            X = [X tau];
            labels = {labels, 'tau'};
        end

    case 'pure_h2'

        % Special set for pure hydrogen
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

    case 'ammonia'

        major_sp = {'o2', 'nh3', 'oh', 'no', 'h2', 'no2', 'o', 'h'};
        id_keep = [];
        count = 1;
        for i = 1 : length(major_sp)
            lab = major_sp{i};
            for j = 1 : length(labels)
                current_label = strrep(labels{j}, ' ', '');
                if strcmp(lab, current_label) == true
                    id_keep(count) = j;
                    count = count + 1;
                end
            end
        end
         
         X = [X_data(:,start) X_data(:,id_keep).^0.5];

end

