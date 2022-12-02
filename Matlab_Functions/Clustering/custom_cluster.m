function [idx, infos] = custom_cluster(k, X, labels, alg, opt)
% This function will cluster the data in k clusters via the selected
% unsupervised learning algorithm 'alg' with the selected options 'opt'
% INPUTS:
%   k = int (number of clusters)
%   X = data matrix (centered and scaled)
%   labels = string vector with labels variables
%   alg = string (algorithm to cluster, available 'k-means', 'lpca', 'gmm'
%   opt = struct (struct variable with the options)

% K-Means options
% opt.Start: 'uniform', 'sample', 'plus', 'center' (see Matlab doc)
% opt.MaxIter: int (max number of iterations)
% K-Means infos
% infos.C = array (cluster centroids)
% infos.sumd = float (sum of squared distances)

% LPCA options
% opt.Precod = bool (preconditioning on mixture fraction)
% opt.f = array (Mixture fraction field if opt.Precond is active)
% opt.StopRule = str (Stopping rule for PCA, available are 'eigs', 'tot_var')
% opt.Inputs = float (Inputs for LPCA, e.g. amount of variance or n eigenvectors, see local_pca_lt.m)

% Check dimensions of X
[np, nv] = size(X);

switch alg
    case 'k-means'

        % For reproducibility
        rng('default');

        % Check for options
        % Initialization
        if isfield(opt, 'Start')
            init = opt.Start;
            opt_avail = {'uniform', 'sample', 'plus', 'cluster'};

            if ~any(strcmpi(init, opt_avail))
                fprintf('Start option not available, will be set to uniform')
                opt.Start = 'uniform';
            end
        else
            fprintf('No initialization options selected for K-Means, uniform set by default');
            opt.Start = 'uniform';
        end

        % Max iterations
        if isfield(opt, 'MaxIter') == true
            itmax = opt.MaxIter;
        else
            fprintf('Max iterations not speciefied, will be set to 1000');
            itmax = opt.MaxIter;
        end

        % K-Means clustering
        % Perform clustering
        [idx, C, sumd] = kmeans(X, k, 'Start', opt.Start, 'MaxIter', opt.MaxIter);

        infos.C = C;
        infos.sumd = sumd;
        infos.title = 'Infos from k-means clustering';

%%%%%%%%%%%%%%%%%%%%%%%%%%%% LPCA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'lpca'

        % Check for preconditioning on mixture fraction
        if isfield(opt, 'Precond') == true
            if opt.Precond == true
                fprintf('Preconditioning on mixture fraction for LPCA \n');
                if isfield(opt, 'f') == true
                    f = opt.f;
                else
                    warning('Mixture fraction not specified in options, standard LPCA will be used');
                    opt.Precond = false;
                end
            else
                fprintf('No preconditioning on mixture fraction');
            end
        else
            fprintf('No options on preconditioning on mixture fraction, raw LPCA will be performed by default');
        end

        % Check for stop rule and options
        if isfield(opt, 'StopRule') == true 
            % Check if stop rule selected is available in the list
            if ~any(strcmp(opt.StopRule, {'var', 'ind_var', 'broken_stick', 'eigs', 'large_eigs'}))
                fprintf('No available stop rule was selected, LPCA with 2 eigenvectors will be run by default');
                opt.StopRule = 'eigs';
                opt.Inputs = 2;

            else
                % Safety checkings on variance
                if strcmp(opt.StopRule, 'var') == true
                    if opt.Inputs > 1
                        warning('Amount of variance should be less than 1. Set to 0.99 by default');
                        opt.Inputs = 0.99;
                    elseif isfield(opt, 'Inputs') == false
                        fprintf('Amount of variance not speciefied. Will be set to 0.99 by default');
                        opt.Inputs = 0.99;
                    end

                    opt.id = 1;

                 % Safety checkings on eigenvectors
                 if strcmp(opt.StopRule, 'eigs') == true
                     if opt.Inputs > nv
                         warning('Number of eigenvectors greater than number of variables. Will be set to 2 by default');
                         opt.Inputs = 2;
                     elseif isfield(opt, 'Inputs') == false
                         fprintf('Number of eigenvectors not specified, will be set to 2 by default');
                         opt.Inputs = 2;
                     end
                     opt.id = 4;
                 end
               end
            end
            
        % If stop rule not specified
        else
            fprintf('Stop rule not specified, will be set to 99 percent of variance by default');
            opt.StopRule = 'var';
            opt.Inputs = 0.99;
            opt.id = 1;
        end

        % Perform FPCA
        if opt.Precond == true
            if opt.Inputs < 1
                fprintf('Only eigenvectors choice is available for FPCA. Will be set to 2 by default');
                opt.Inputs = 2;
            end
            idx = local_pca_lt_2(X, f, labels, k, opt.Inputs, opt);

        else

            [idx, infos.eigvec, infos.eps_rec, infos.nz_X_k, infos.U_scores, infos.rec_data, infos.gamma_pre, infos.C, infos.nz_idx_clust, infos.X_ave_pre, infos.min_var_expl_clust, infos.sq_rec_err, infos.var_exp, infos.X_ave_clust] = ...
                local_pca_lt(X, labels, k, nv-2, [], 0, 0, opt.StopRule, opt.Inputs);

        end

end