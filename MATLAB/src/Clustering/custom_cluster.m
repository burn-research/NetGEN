function [idx, infos] = custom_cluster(k, X, labels, alg, opt)
%
% function [idx, infos] = custom_cluster(k, X, labels, alg, opt)
%
%
% This function will cluster the data in k clusters via the selected
% unsupervised learning algorithm 'alg' with the selected options 'opt'
% INPUTS:
%   k = int (number of clusters)
%   X = data matrix (centered and scaled)
%   labels = string vector with labels variables
%   alg = string (algorithm to cluster, available 'k-means', 'vqpca')
%   opt = struct (struct variable with the options)
%
% OUTPUTS
%
%   idx   = cluster labels vector
%   infos = struct variable with clustering outputs depending on algorithm
%           used
%
%%%%%%%%% Available options for kmeans %%%%%%%%%%
%
%   opt.Start   = string {'uniform', 'sample', 'plus', 'cluster'}
% (initialization method)
%
%   opt.MaxIter = int (maximum number of iterations)
%
%%%%%%%%% Available options for Local PCA %%%%%%%%
%   
%   opt.Precond = bool (if true, data will be preconditioned on mixture
%   fraction f)
%
%   opt.StopRule = string {'var', 'ind_var', 'broken_stick', 'eigs', 'large_eigs'}
%   (defines stopping rule used for VQPCA)
%
%       'var'          = amount of variance retained (float 0 < var < 1)
%       'ind_var'      = individual variance that must be retained at least for
%                        all the variables (float 0 < ind_var < 1)
%       'broken_stick' = see pca_lt.m
%       'eigs'         = number of eigenvectors to retain
%       'large_eigs'   = check largest eigenvectors with hard threshold
%
%   opt.Inputs   = check vqpca.m

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%% VQPCA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'vqpca'

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
            idx = fpca(X, f, labels, k, opt.Inputs, opt);

        else

            [idx, infos] = vqpca(X, k, opt.StopRule, opt.Inputs, opt);

        end

end