% LOCAL PCA - Performs local pca on a data set.
% 
% function [idx, infos] = local_pca_lt(X, k, stop_rule, input, opt)
% INPUTS
% 
% X             = Matrix of state variables
% k             = Number of clusters (integer)
% stop_rule     = Refer to pca_lt.m (integer)
% input         = Depending on the stop rule (e.g. variance fraction)
% opt           = Struct variable containing different options (see
%                 documentation below)
%
% OPTIONAL:
%
%   opt = Structure variable of matlab containing several options as
%   fields. As an example opt.Scaling = 'auto'
%
%   Available options are:
%   
%   opt.Scaling = (string) {'auto', 'range', 'pareto', 'vast', 'max', 'level'} see
%   scale.m for more info (auto by default)
%
%   opt.Center  = (integer) Select 1 for centering and 0 for no centering
%   (1 by default)
%
%   opt.MaxIter = (integer) Maximum number of iterations (600 by default)
%   
%   opt.EpsRecMin = (float) Minimum threshold of ratio between cluster centroids
%   across two successive iterations (1.0e-6 default)
%
%   opt.Algorithm = (string) Available are VQPCA and FPCA. In the second
%   case more fields are required (VQPCA by default)
%
%   opt.Init = (string) Initialization method. Available are:
%               {'random', 'uniform1', 'PCA', 'bestDB', 'uniform2',
%               'uniform3'}. Refer to initialize_centroids.m for more info
%               (uniform1 by default)
%
%   opt.C = (array) Array of initial centroids for custom initialization
%   (empty by default)
%
%   opt.InitEigs = (integer) 1 for standard (from ones and identity
%   matrices) 2 for an initial PCA (1 by default)
%
%   opt.EigStart = (integer) If opt.InitEigs = 2, you must specify the
%   number of eigenvectors for the first iteration. Must be lower than
%   number of columns (2 by default)
%
%   % ONLY IN CASE FPCA IS CHOSEN OTHER FIELDS REQUIRED ARE:
%
%   opt.F = (vector) Array containing the point-wise mixture fraction
%
%   opt.FStoich = (float) Stoichiometric mixture fraction
%
%   OUTPUTS: 
%
%       idx = vector of clustering labels
%
%       infos = struct array with the fields listed below
%
% infos.C = C; 
% infos.eigenvectors = eigvec;  
% infos.eigenvalues  = eigval; 
% infos.ReconstructionError = eps_rec_new; 
% infos.MinimumReconstructionError = eps_rec_minimum;
% infos.MaximumReconstructionError = eps_rec_maximum;
% infos.WeightedReconstructionError = eps_rec_weighted;
% infos.ExplainedVariance = var_exp;
% infos.LocalDimensionality = n_eigs;
% infos.UScores = U_scores;
% infos.CPUTime = cputime;
% infos.WScores = W_scores;
% infos.Gamma = gamma;
% infos.X_ave_clust = X_ave_clust;
% infos.Loadings = loadings;
%
%
% This routine peforms local PCA on a data set. Two different partitions
% algorithms can be selected: a supervised (FPCA) and unsupervised (VQPCA)
% one. 
% 1) When using FPCA data are partitioned into bins of mixture fraction
% and then PCA is performed locally in each one of the bin. 
% 2) When VQPCA is performed, the data are assigned into clusters based on
% their reconstruction distance, i.e. the diference between the
% actual p-dimesnional point and its q-dimesnional approximation. The
% reconstruction error is then defined as eGSRE=X-Xq. The procedure is then
% iterative and it follows the following steps: 
% (a) Initialization: the cluster centroids are randomly chosen from the
% data set and the covariance matrix is initialized to the identity matrix
% for each cluster
% b) Partition: each observation from the sample is assigned to a cluster
% using the squared reconstruction distance given by eGSRE
% c) Update: the clusters??? centroids are updated on the basis of
% partitioning
% d) Local PCA: PCA is performed in each disjoint region of the sample.
% e) Steps 2-4 are iterated until convergence is reached.
% The goodness of reconstruction given by VQPCA is measured with respect to
% the mean variance in the data as eGSRE,n=E(eGSRE)/E(var(X)). If the auto
% scaling criterion is adopted the mean variance of X after scaling equals
% one and the reconstruction error metric becomes eGSRE,n=E(eGSRE).
% Convergence can be judged using the following criteria:
% a) The normalized global scaled reconstruction error, eGSRE,n, is below a
% specific threshold, eGSRE,n.
% b) The relative change in clusters??? centroids between two successive
% iterations is below a fixed threshold, i.e. 10-8.
% c) The relative change in eGSRE,n between two successive iterations is
% below a fixed threshold, i.e. 10-8. Requirements b) and c) are
% particularly useful if an explanatory analysis on the performances of
% VQPCA in terms of eGSRE,n is of interest. In this case, requirement a)
% can be relaxed and the variation of eGSRE,n as a function of the number
% of eigenvalues and clusters can be analyzed by enforcing requirements b)
% and c). Otherwise, all the three conditions can be used and an iterative
% procedure for the determination of the number of eigenvalues required to
% achieve a fixed eGSRE,n could be employed. Staring with q=1, the number
% of eigenvalues can be increased progressively until the desired error
% level is reached.
% Main function

function [idx, infos] = ...
    vqpca(X, k, stop_rule, input, opt)

% Number of rows and columns
[rows, columns] = size(X);

% Set font size
set(0,'defaultaxesfontsize',20);
set(0,'defaulttextfontsize',20);

%% Convergence parameters

% Maximum number of iterations
if isfield(opt, 'MaxIter') == false
    iter_max = 200;
else
    iter_max = opt.MaxIter;
    if isinteger(iter_max) == false
        error('Max number of iterations must be an integer');
    end
end

% Eps rec threshold
if isfield(opt, 'EpsRecMin') == false
    eps_rec_min = 1.0e-6;
else
    eps_rec_min = opt.EpsRecMin;
end

% Convergence indicators initialization
convergence = 0;
iter = 0;
eps_rec = 1.0;
a_tol = 1.0e-10;
r_tol = 1.0e-10;

% Initial time
t0 = cputime;

%% CENTER AND SCALE THE DATA
% Select centering and scaling criteria
% Before clustering we do not center, because  we remove the centroid
% inside the clustering we are centering (we do not center twice)

% Pre-centering option
if isfield(opt, 'Center') == false
    pre_cent_crit = 1;
    disp('Centering option not specified, data will be centered by default');
else
    pre_cent_crit = opt.Center;
    if isinteger(pre_cent_crit) == false
        warning('Specify 1 for centering, 0 for no centering. Data will be centered');
        pre_cent_crit = 1;
    end
end

% Pre-scaling option
if isfield(opt, 'Scaling') == false
    pre_scal_crit = 1;
    opt_3 = 'auto';
    disp('Data will be scaled using the auto scaling by default');
else
    scaling_list = {'auto', 'range', 'pareto', 'vast', 'level', 'max', 'no'};
    if ismember(opt.Scaling, scaling_list) == false
        warning('Invalid scaling criterion specified. Data will be scaled with auto scaling');
    else
        if strcmp(opt.Scaling, 'auto')
            pre_scal_crit = 1;
        elseif strcmp(opt.Scaling, 'range')
            pre_scal_crit = 2;
        elseif strcmp(opt.Scaling, 'pareto')
            pre_scal_crit = 3;
        elseif strcmp(opt.Scaling, 'vast')
            pre_scal_crit = 4;
        elseif strcmp(opt.Scaling, 'level')
            pre_scal_crit = 5;
        elseif strcmp(opt.Scaling, 'max')
            pre_scal_crit = 6;
        elseif strcmp(opt.Scaling, 'no')
            pre_scal_crit = 0;
        end
    end

    opt_3 = opt.Scaling;
end

% Scale and center
[cent_X, X_ave] = center(X, pre_cent_crit);
[scal_X, X_gamma] = scale(cent_X, X, pre_scal_crit);
gamma_pre = X_gamma;
X_ave_pre = X_ave;

% Save scaling factors in infos
infos.gamma_pre = gamma_pre;
infos.X_ave_pre = X_ave_pre;

% Choose the clustering criterion
if isfield(opt, 'Algorithm') == false
    VQ = 1;
    disp('Unsupervised VQPCA will be applied');
else
    alg = opt.Algorithm;
    if strcmp(alg, 'VQPCA')
        VQ = 1;
    elseif strcmp(alg, 'FPCA')
        VQ = 2;
    end

    if VQ ~=1 || VQ ~= 2
        warning('Invalid choice of algorithm, VQPCA will be selected');
    end
end

if VQ == 1
    VQ_name = 'VQPCA';
else
    VQ_name = 'FPCA';
end

%% 1) INITIALIZATION

% Centroids initialization
if VQ == 1
    if isfield(opt, 'C') == false
        disp('Centroids not specified. Initialization required');
        if isfield(opt, 'Init') == false
            init = 1;
            disp('Initialization method not specified. Uniform1 will be chosen');
            opt.Init = 'uniform1';
        else
            init = opt.Init;
        end

        % Check function initialize_centroids.m for more info
        C = initialize_centroids(scal_X, k, opt);
    else
        C = opt.C;
        disp('Centroids were already specified');
        [nn,mm] = size(C);
        if nn~= k || mm~=columns
            error('Centroids have wrong dimensions');
        end
    end
    
elseif VQ == 2
    opt_3 = 'Data partitioned on F';
    %The stoichiometric Z should be known.
end

%% Initialization of eigenvectors
n_eigs_max = columns-1;
if isfield(opt, 'InitEigs') == false
    init_eigs = 1;
else
    init_eigs = opt.InitEigs;
    if init_eigs ~= 1 || init_eigs ~= 2
        warning('InitEigs specified is not valid. Eigenvectors will be initialized as ones vector');
    end
end

% Number of eigenvectors for start
if isfield(opt, 'EigsStart') == false
    n_eigs = 2;
else
    n_eigs = opt.EigsStart;
    if n_eigs > n_eigs_max
        warning('Too many eigenvectors!');
        n_eigs = n_eigs_max;
    end
end

% Initialization as ones and zeros vectors
if init_eigs == 1
    fprintf('\nInitialization of the eigenvectors as ones and zeros vectors \n');
    eigvec = cell(k, 1);
    gamma = cell(k, 1);
    for j = 1 : k
        eigvec{j} = eye(columns, n_eigs);
        gamma{j} = ones(1, columns);
    end
    
% Initialization from previous PCA
elseif init_eigs == 2
    fprintf('\n Initialization of the eigenvectors from preliminary PCA \n');
    eigvec = cell(k, 1);
    gamma = cell(k, 1);
    for j = 1 : k
        [sort_eigval, sort_eigvec, ret_eigval, ret_eigvec, n_eig, U_scores, W_scores, gamma{j}] ...
            = pca_lt(X, pre_cent_crit, pre_scal_crit, stop_rule, input);
        eigvec{j} = ret_eigvec;
    end

% Initialization from the starting clustering, valid only if Peppe init is
% used
elseif init_eigs == 3
    fprintf('\n Initialization of the eigenvectors from first PCA \n');
    eigvec = cell(k, 1);
    gamma = cell(k, 1);
    X_clust = clustering(scal_X, idx);
    for j = 1 : k
        [sort_eigval, sort_eigvec, ret_eigval, ret_eigvec, n_eig, U_scores, W_scores, gamma{j}] ...
            = pca_lt(X_clust{j}, 0, 0, stop_rule, input);
        eigvec{j} = ret_eigvec;
    end
end

% CREATE A DIRECTORY FOR THE OUTPUTS
thetime = clock;
thetimestr = datestr(thetime);
dirname = [VQ_name '_n_clust_' num2str(k) '_n_eig_' num2str(n_eigs) '_' ...
    thetimestr(1:11) '_' num2str(thetime(4)) '_' num2str(thetime(5)) '_' ...
    num2str(floor(thetime(6)))];
mkdir(dirname)
cd(dirname)
fid = fopen('output.out','w');
fprintf(fid, 'Output file \n');
fprintf(fid, '\nClusters centroids initialization: %s', opt_3);
fprintf(fid, '\n');
fprintf(fid, '\nInitial clusters centroids \n');
if VQ == 1
    for j = 1 : k 
        for l = 1 : columns
            fprintf(fid,'%d \t', C(j, l));
        end
        fprintf(fid,'\n');
    end
end

%% 2) PARTITION

% Partition can be performed in two ways:
% 1) Iteratively, using the unsupervised quantization algorithm
% 2) By conditioning the data using a supervised quantization algorithm,
% i.e. by conditioning the data on the mixture fraction

% Select centering and scaling criteria for PCA
cent_crit = 1;
scal_crit = 0;

min_var_expl_clust = ones(k,1);

while ((convergence == 0) && (iter < iter_max))
    
    C_convergence = 0;
    eps_rec_convergence = 0;   
    fprintf('\nIteration n. %d, convergence %d \n', iter, convergence);  

    sq_rec_err = zeros(rows, k);
    sq_rec_err_pen = zeros(rows, k);
    
    if ((VQ==2 && iter == 0))

        % Check if opt.F exists
        if isfield(opt, 'F') == false
            error('FPCA was selected but no mixture fraction was specified. Please specify mixture fraction as opt.F = array');
        else
            F = opt.F;
        end

        % Check if stoichiometric mixture fraction was specified
        if isfield(opt, 'FStoich') == false
            error('opt.FStoich not specified. Please specify it as float');
        else
            F_stoich = opt.FStoich;
        end
        F_min = min(F);
        F_max = max(F);

        %F_stoich = 0.4375;% DNS CO/H2 Flame
        %F_stoich = 0.351; % Flame F

        [nz_X_k, nz_idx_clust] = condition(scal_X, F, k, F_min, F_max, F_stoich);
        C = zeros(k, columns);
        for j = 1 : k
            C(j, :) = mean(nz_X_k{j}, 1);
            [sort_eigval{j}, sort_eigvec{j}, eigval{j}, eigvec{j}, n_eig{j}, ...
                U_scores{j}, W_scores{j}, gamma{j}, scaled_data{j}, ...
                rec_data{j}] = pca_lt(nz_X_k{j}, cent_crit, scal_crit, stop_rule, input);
            for l = 1 : columns
                fprintf(fid,'%d \t', C(j, l));
            end
            fprintf(fid,'\n');
        end

    end

    % Evaluate reconstruction error
    if isfield(opt, 'CustomError')
        sq_rec_err = custom_rec_err(scal_X, C, gamma, eigvec, opt);
    else
        % Evaluate the squared mean reconstruction error
        for j = 1 : k
            D = diag(gamma{j});
            C_mat = repmat(C(j, :), rows, 1);
            rec_err_os = (scal_X - C_mat - (scal_X - C_mat) * D^-1 * eigvec{j} * eigvec{j}' * D);
            sq_rec_err(:, j) = sum(rec_err_os.^2, 2);
        end
    end
    
    % Classical reconstruction error
    [rec_err_min, idx] = min(sq_rec_err, [], 2);
    rec_err_min_rel = (rec_err_min);
    
    % Squared reconstruction error of each observation, in each cluster
    [sq_err_clust, nz_idx_clust, k] = partition(rec_err_min, idx, k);
    
    % Evaluate the global mean error
    eps_rec_new = mean(rec_err_min_rel);
      
    if VQ == 1
        
        % Partition the data into clusters
        [nz_X_k, nz_idx_clust, k]  = partition(scal_X, idx, k);
        
        if length(nz_X_k) ~= length(nz_idx_clust)
            warning('There is a problem');
            pause;
        end
    end
    
    if isempty(nz_X_k) == true
        warning('Empty partition of data, using other function to partition');
        nz_X_k = clustering(scal_X, idx);
        nz_idx_clust = clustering([1:1:length(scal_X)]', idx);
    end
    
    % Evaluate the relative recontruction errors in each cluster
    rec_err_min_rel_k = cell(k, 1);
    for j = 1 : k
        rec_err_min_rel_k{j} = rec_err_min_rel(nz_idx_clust{j}, 1);
    end
    
    % Evaluate the mean error in each cluster
    eps_rec_new_clust = zeros(k, 1);
    size_clust = zeros(k, 1);
    rec_err_clust = zeros(k, 1);
    wi = 0;
    
    for j = 1 : k
        eps_rec_new_clust(j) = mean(rec_err_min_rel_k{j});
        size_clust(j) = size(nz_idx_clust{j}, 1);
        rec_err_clust(j) = sum(sq_err_clust{j})/size_clust(j);   % NMRSE
        wi = eps_rec_new_clust(j)*size_clust(j) + wi;
    end
    
    % Evaluate different reconstruction errors
    eps_rec_minimum = min(eps_rec_new_clust);
    eps_rec_maximum = max(eps_rec_new_clust);   
    eps_rec_weighted = wi/rows;

    fprintf('\nGlobal mean recontruction error at iteration n. %d equal to %d \n', iter, eps_rec_new);
    fprintf('\nLocal mean recontruction error at iteration n. %d \n', iter);
    for j = 1 : k
        fprintf('%d \t', eps_rec_new_clust(j));
    end
    fprintf('\n');
    
    if VQ == 1
        
        % 3) EVALUATE NEW CLUSTERS' CENTROIDS
    
        C_new = zeros(k, columns);        
        for j = 1 : k
            C_new(j, :) = mean(nz_X_k{j}, 1);
        end
        
        eps_rec_var = abs((eps_rec_new  - eps_rec) / eps_rec_new);
        fprintf('\nReconstruction error variance equal to %d \n', eps_rec_var);

        if ((eps_rec_var < r_tol) && (eps_rec_new > eps_rec_min) ...
                && (n_eigs < n_eigs_max)) 
            n_eigs = n_eigs + 1;
            fprintf('\n Cluster %d dimension increased to %d \n', j,  n_eigs);
        end

        % Judge convergence: clusters centroid and relative reconstruction
        % error
    
        if (eps_rec_var < r_tol)
            eps_rec_convergence = 1;
        end
        if (size(C) == size(C_new))
            C_var = abs((C_new - C) ./ (C_new + a_tol));
            if (C_var(:, :) < r_tol)
                C_convergence = 1;
            end
        end
        if ((iter > 1) && (C_convergence == 1) && (eps_rec_convergence == 1))
            convergence = 1;
            fprintf('\nConvergence reached in %d iterations \n', iter);
        end

        % Update recontruction error and cluster centroids
        C = C_new;
        eps_rec = eps_rec_new;
        
        % 4) PERFORM LOCAL PCA
   
        % Initialization of cell arrays
        sort_eigval = cell(k, 1);
        sort_eigvec = cell(k, 1);
        eigval = cell(k, 1);
        eigvec = cell(k, 1);
        n_eig = cell(k, 1);
        U_scores = cell(k, 1);
        W_scores = cell(k, 1);
        gamma = cell(k, 1);
        scaled_data = cell(k, 1);
        rec_data = cell(k, 1);
        X_ave_clust = cell(k,1);
        
        % Minimum individual variance explained in each cluster
        min_var_expl_clust = zeros(k,1);
        
        % Variance of the single variables explained in each cluster
        var_exp = cell(k,1);
        
        % Perform PCA in each cluster
        for j = 1 : k
                        
            [sort_eigval{j}, sort_eigvec{j}, eigval{j}, eigvec{j}, n_eig{j}, ...
                U_scores{j}, W_scores{j}, gamma{j}, scaled_data{j}, rec_data{j}, X_ave_clust{j}] ...
                = pca_lt(nz_X_k{j}, cent_crit, scal_crit, stop_rule, input);
            
            
            % Calculate the explained variance for each data
            loadings = zeros(columns, length(eigval{j}));
            cov_data = cov(nz_X_k{j});
            for i = 1 : length(eigval{j})
                for l = 1 : columns
                    loadings(l, i) = (eigvec{j}(l,i)*eigval{j}(i)^0.5)/(cov_data(l,l)^0.5);
                end
            end
            
            var_expl_x = zeros(columns, 1);
            for i = 1 : columns
                var_expl_x(i) = 0;
                for l = 1 : length(eigval{j})
                    var_expl_x(i) = var_expl_x(i) + loadings(i,l)^2;
                end
            end
            
            [min_var_expl_clust(j), ~] = min(var_expl_x);
            var_exp{j} = var_expl_x;
                        
        end
        
        iter = iter + 1;
        
    elseif VQ == 2
        convergence = 1;    
    end
end

if (convergence == 0)
    fprintf('\nConvergence not reached in %d iterations \n', iter);
end

% MIXTURE FRACTION PARTITION
if VQ == 2
    F_k = cell(k, 1);
    unscaled_nz_X_k = cell(k, 1);
    uncentered_nz_X_k = cell(k, 1);
    for j = 1 : k
        F_k{j} = F(nz_idx_clust{j}, 1);
        unscaled_nz_X_k{j} = unscale(nz_X_k{j}, X_gamma);
        uncentered_nz_X_k{j} = uncenter(unscaled_nz_X_k{j}, X_ave(nz_idx_clust{j}, :));    
    end
end

% ORIGINAL DATA AND RECOVERED DATA RECONSTRUCTION

rec_scal_X = zeros(rows, columns);
rec_scal_X_hat = zeros(rows, columns);
for j = 1 : k
    rec_scal_X(nz_idx_clust{j}, :) = nz_X_k{j};
    rec_scal_X_hat(nz_idx_clust{j}, :) = rec_data{j};
end
rec_cent_X_data = unscale(rec_scal_X_hat, X_gamma);
rec_X_data = uncenter(rec_cent_X_data, X_ave);

% Check that the reconstructed original data are the same of the original
% data

if (abs(scal_X(:, :) - rec_scal_X(:, :)) > a_tol)
    error('The reconstructed data are non equal to the original data');
end

% CPU TIME
overall_cpu_time = cputime - t0;

% WRITE THE OUTPUT FILE

fprintf(fid, '\nTotal number of clusters equal to %d \n', k);
fprintf(fid, '\nTotal number of iterations equal to %d \n', iter);
fprintf(fid, '\nRelative recontruction error equal to %d \n', eps_rec_new);
fprintf(fid, '\nRelative recontruction errors in each cluster \n');
for j = 1 : k
    fprintf(fid, '%d \t', eps_rec_new_clust(j)); 
end
fprintf(fid, '\n');
fprintf(fid, '\nNumber of eigenvalues in each cluster \n');
for j = 1 : k
    fprintf(fid, '%d \t', n_eigs); 
end
fprintf(fid, '\n');
fprintf(fid, '\nTotal CPU time equal to %d s \n', overall_cpu_time);
% for j = 1 : k
%     fprintf(fid, '\nCluster n. %d size = ', j);
%     fprintf(fid, '%d \n', size_clust(j));
% end
fprintf(fid, '\nFinal clusters centroids \n');
for j = 1 : k 
    for l = 1 : columns
    fprintf(fid,'%d \t', C(j, l));
    end
    fprintf(fid,'\n');
end

% SAVE DATA 

save nz_X_k nz_X_k
% save F_k F_k
save rec_X_data rec_X_data
save idx.out idx -ASCII -DOUBLE
save C.out C -ASCII -DOUBLE
save eigval eigval
save eigvec eigvec
save U_scores U_scores

%% Save infos to the infos variable
infos.C = C;
infos.eigenvectors = eigvec;
infos.eigenvalues  = eigval;
infos.ReconstructionError = eps_rec_new;
infos.MinimumReconstructionError = eps_rec_minimum;
infos.MaximumReconstructionError = eps_rec_maximum;
infos.WeightedReconstructionError = eps_rec_weighted;
infos.ExplainedVariance = var_exp;
infos.LocalDimensionality = n_eigs;
infos.UScores = U_scores;
infos.CPUTime = cputime;
infos.WScores = W_scores;
infos.Gamma = gamma;
infos.X_ave_clust = X_ave_clust;
infos.Loadings = loadings;
infos.RecData = rec_data;
infos.NzIdxClust = nz_idx_clust;

cd ..
                
    