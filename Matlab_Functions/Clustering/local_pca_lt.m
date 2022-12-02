% LOCAL PCA - Performs local pca on a data set.
% 
% function [idx, eigvec, eps_rec, nz_X_k, U_scores, rec_data, gamma_pre, C, nz_idx_clust, X_ave_pre, min_var_expl_clust, sq_rec_err, var_exp, X_ave_clust] = local_pca_lt(X, names, k, n_eigs, C, pre_cent_crit, pre_scal_crit, stop_rule, input)
% INPUTS
% 
% X             = Matrix of state variables
% F             = Mixture fraction
% names         = State variables names (MUST be cellstr format)
% k             = Number of clusters (integer)
% n_eigs        = Max number of eigenvectors (integer)
% pre_cent_crit = 1 for centering, 0 for no centering (integer)
% pre_scal_crit = Refer to scale.m (integer)
% stop_rule     = Refer to pca_lt.m (integer)
% input         = Depending on the stop rule (e.g. variance fraction)
%
% OPTIONAL:
% C            = cluster centroids initialization
%
% OUTPUTS
%
% idx               = Vector of indexes denoting the cluster number for
%                     each of the original data points
% eigvec            = cell k x n_eigs containing eigenvectors in each
%                     clusters
% eps_rec           = average reconstruction error in each cluster, vector
%                     k x 1
% nz_X_k            = Cell matrix of the partitioned data, i.e. each  cell
%                     is a cluster of points
% U_scores          = Cell k x 1 containing the U_scores in each cluster
% rec_data          = reconstructed data
% gamma_pre         = scaling factor for each cluster and each variable
% C                 = centroids
% nz_idx_clust      = cell k x 1 containing the label of each point
% X_ave_pre         = pre-centerd and pre-scaled data matrix
% min_var_expl_clust = minimum variance explained in each cluster (refers
%                      to the variance explained for the single variables)
% sq_rec_err        = n_rows x k matrix, element (i,j) is the
%                     reconstruction error for point i in cluster j (refer to db_PCA.m)
% var_exp           = variance explained for each variable in each cluster
% X_ave_clust       = unscaled and uncentered data in each cluster
% 
% 
% F_k               = Cell matrix of mixture storing the values of mixture
%                     fraction in each bin identified by LPCA
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

function [idx, eigvec, eps_rec, nz_X_k, U_scores, rec_data, gamma_pre, C, nz_idx_clust, X_ave_pre, min_var_expl_clust, sq_rec_err, var_exp, X_ave_clust] = ...
    local_pca_lt(X, names, k, n_eigs, C, pre_cent_crit, pre_scal_crit, stop_rule, input)

[rows, columns] = size(X);

% INPUTS
% Choose the number od clusters
plots = 0;

% Set font size
set(0,'defaultaxesfontsize',16);
set(0,'defaulttextfontsize',18);
set(0, 'defaulttextfontweight', 'Bold')

% Convergence indicators initialization
convergence = 0;
iter = 0;
iter_max = 200;
eps_rec = 1.0;
eps_rec_min = 1.0e-06;
a_tol = 1.0e-10;
r_tol = 1.0e-10;

% Initial time
t0 = cputime;

% CENTER AND SCALE THE DATA
% Select centering and scaling criteria
% Before clustering we do not center, because  we remove the centroid
% inside the clustering we are centering (we do not center twice)

[cent_X, X_ave] = center(X, pre_cent_crit);
[scal_X, X_gamma] = scale(cent_X, X, pre_scal_crit);
gamma_pre = X_gamma;
X_ave_pre = X_ave;

% Choose the clustering criterion
VQ = 1;
if VQ == 1
    VQ_name = 'VQPCA';
else
    VQ_name = 'FPCA';
end

% 1) INITIALIZATION

% Centroids initialization

if VQ == 1
    
    % If no centroids initialization is provided specify the default
    % initialization method with the variable init
    % init = 1: uniform
    % init = 2: random
    % init = 3: from FPCA
    % init = 4: from global PCA sampling on U scores
    % init = 5: best from random choice of 10 iterations
    % init = 6: peppe uniform
    % init = 7: sort on temperature
    
    if isempty(C) == true
        fprintf('\nClusters centroids not specified, initialization required \n');
        init = 6;
        if init == 1
            C_int = linspace(1, rows, k+2);
            C = scal_X(round(C_int(2:k+1)), :);
            opt_3 = 'uniform';
        elseif init == 2
            C_int = randomsample(rows, k);
            C = scal_X(C_int, :);
            opt_3 = 'random';
        elseif init == 3
            opt_3 = 'from FPCA';
        elseif init == 4
            
            % Random initialization from PC scores
            [~, ~, ~, ~, n_eig, U_scores, ~, ~, ~, ~] = ...
                pca_lt(X, pre_cent_crit, pre_scal_crit, 1, 0.95);
            
            % Range of the U_score matrix
            inter = linspace(min(U_scores(:,1)), max(U_scores(:,1)), k+1);
            
            % Pick all the values in between the intervals and randomly
            % select an index from them
            C_int = zeros(k,1);
            for i = 1 : length(inter) - 1
                samples = find(U_scores(:,1) > inter(i) & U_scores(:,1) < inter(i+1));
                ind = randsample(length(samples), 1);
                C_int(i) = ind;
            end
            
            % Centroids initialization
            C = scal_X(C_int, :);
            opt_3 = 'U score range';
        
        % Initialize with best random of 10 initial iterations
        elseif init == 5
            it_init = 100;
            fprintf('Initializing from the best of initial %d iterations /n', it_init);
            rec_err_init = zeros(it_init, 1);
            C_init = cell(it_init, 1);
            db_init = zeros(it_init, 1);
            for i = 1 : it_init
                C_int = randomsample(rows, k);
                C = scal_X(C_int, :);
                C_init{i} = C;
                
                % Init eigenvectors
                eigvec = cell(k, 1);
                gamma = cell(k, 1);
                sq_rec_err = zeros(rows, k);
                for j = 1 : k
                    eigvec{j} = eye(columns, n_eigs);
                    gamma{j} = ones(1, columns);
                
                    % Evaluate reconstruction error at starting point
                    D = diag(gamma{j});
                    C_mat = repmat(C(j, :), rows, 1);
                    rec_err_os = (scal_X - C_mat - (scal_X - C_mat) * D^-1 * eigvec{j} * eigvec{j}' * D);
                    sq_rec_err(:, j) = sum(rec_err_os.^2, 2);        
                end
                
                
                rec_err_init(i) = mean(min(sq_rec_err, [], 2));
                [~, idx_init] = min(sq_rec_err, [], 2);
                
                db_init(i) = db_pca_mod(X, idx_init, stop_rule, input);
                
            end
            
            % Pick minimum from all the iterations
            [~, id_min] = min(rec_err_init);
            C = C_init{id_min};
            opt_3 = 'Best_random';
            
        elseif init == 6
            
            spacing = ceil(rows/k);
            idx = zeros(rows,1);
            
            for i = 1 : rows
                idx(i) = ceil(i/spacing);
            end
            
            data_part = clustering(scal_X, idx);
            
            C = zeros(k,columns);
            for i = 1 : length(data_part)
                C(i,:) = mean(data_part{i}, 1);
            end
            
            opt_3 = 'peppe';

        elseif init == 7

            [~,index] = sort(X(:,1), 'ascend');
            scal_X = scal_X(index,:);
            C_int = linspace(1, rows, k+2);
            C = scal_X(round(C_int(2:k+1)), :);
            opt_3 = 'uniform';
            
        end
    else
        fprintf('\nClusters centroids initialized \n');
        opt_3 = 'Initialized';
    end
elseif VQ == 2
    opt_3 = 'Data partitioned on F';
    %The stoichiometric Z should be known.
end

%% Initialization of eigenvectors
init_eigs = 1;
n_eigs_max = columns-1;

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
if (VQ == 1 && init ~= 3)
    for j = 1 : k 
        for l = 1 : columns
            fprintf(fid,'%d \t', C(j, l));
        end
        fprintf(fid,'\n');
    end
end

% 2) PARTITION

% Partition can be performed in two ways:
% 1) Iteratively, using the unsupervised quantization algorithm
% 2) By conditioning the data using a supervised quantization algorithm,
% i.e. by conditioning the data on the mixture fraction

% Select centering and scaling criteria for PCA
cent_crit = 1;
scal_crit = 0;

% Decide to apply the penalty
% 0 = No penalty
% 1 = Variance
% 2 = kurtosis
% 3 = Skewness
% 4 = Kernel

penalty = 0;
min_var_expl_clust = ones(k,1);

while ((convergence == 0) && (iter < iter_max))
    
    C_convergence = 0;
    eps_rec_convergence = 0;   
    fprintf('\nIteration n. %d, convergence %d \n', iter, convergence);  

    sq_rec_err = zeros(rows, k);
    sq_rec_err_pen = zeros(rows, k);
    
    if ((VQ == 2 || init == 3) && iter == 0) 
        F_min = min(F);
        F_max = max(F);
        F_stoich = input('\nInput the value of the stoichiometric mixture fraction: Fstoich = \n');
        %F_stoich = 0.4375;% DNS CO/H2 Flame
        %F_stoich = 0.351; % Flame F
        [nz_X_k, nz_idx_clust] = condition(scal_X, F, k, F_min, F_max, F_stoich);
        C = zeros(k, columns);
        for j = 1 : k
            C(j, :) = mean(nz_X_k{j}, 1);
            [sort_eigval{j}, sort_eigvec{j}, eigval{j}, eigvec{j}, n_eig{j}, ...
                U_scores{j}, W_scores{j}, gamma{j}, scaled_data{j}, ...
                rec_data{j}] = pca_lt(nz_X_k{j}, cent_crit, scal_crit, 1, 0.9);
            for l = 1 : columns
                fprintf(fid,'%d \t', C(j, l));
            end
            fprintf(fid,'\n');
        end

    end

    % Evaluate the squared mean reconstruction error
    for j = 1 : k
        D = diag(gamma{j});
        C_mat = repmat(C(j, :), rows, 1);
        rec_err_os = (scal_X - C_mat - (scal_X - C_mat) * D^-1 * eigvec{j} * eigvec{j}' * D);
        sq_rec_err(:, j) = sum(rec_err_os.^2, 2);
        
%         if penalty ~= 4
%             rec_err_os = (scal_X - C_mat - (scal_X - C_mat) * D^-1 * eigvec{j} * eigvec{j}' * D);
%             sq_rec_err(:, j) = sum(rec_err_os.^2, 2);
%             
%         elseif penalty == 4
%             rec_err_os = (scal_X - C_mat - (scal_X - C_mat) * D^-1 * eigvec{j} * eigvec{j}' * D);
%             [~,n] = size(rec_err_os);
%             sigma = 1000;
%             sq_rec_err(:, j) = (1/(((2*pi)^0.5)*sigma)^n) * exp(-sum(rec_err_os.^2, 2)/(2*sigma^2));
%         end
%         
%         % Variance Penalty
%         if penalty == 1
%             dev_X = scal_X - C_mat;
%             var_X = var(dev_X,0,2);
%             sq_rec_err_pen(:, j) = sq_rec_err(:, j).*var_X;
%         
%         % Kurtosis Penalty
%         elseif penalty == 2
%             dev_X = scal_X - C_mat;
%             kur = kurtosis(dev_X, 0, 2);
%             sq_rec_err_pen(:, j) = sq_rec_err(:, j).*kur;
%         
%         % Skewness Penalty
%         elseif penalty == 3
%             dev_X = scal_X - C_mat;
%             sk = skewness(dev_X, 0, 2);
%             sq_rec_err_pen(:, j) = sq_rec_err(:, j).*sk;
%                         
%         % No penalty
%         else
%             sq_rec_err_pen(:, j) = sq_rec_err(:, j);
%         end      
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

% F_k = cell(k, 1);
% unscaled_nz_X_k = cell(k, 1);
% uncentered_nz_X_k = cell(k, 1);
% for j = 1 : k
%     F_k{j} = F(nz_idx_clust{j}, 1);
%     unscaled_nz_X_k{j} = unscale(nz_X_k{j}, X_gamma);
%     uncentered_nz_X_k{j} = uncenter(unscaled_nz_X_k{j}, X_ave(nz_idx_clust{j}, :));    
% end

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

% PLOTS

% PLOT OF THE ORIGINAL DATA VS RECONTRUCTED DATA

plots = 0;
if plots == 1
    n_clust = num2str(k);
    opt_1 = char(n_clust);
    n_eigs_clust = num2str(n_eigs);
    opt_2 = char(n_eigs_clust);
    for l = 1 : columns
        figure;
        plot(X(:, l), rec_X_data(:, l), 'b*');
        hold on;
        plot(X(:, l), X(:, l), 'r--');
        xlabel(['Observed ', char(names(l))], 'FontWeight', 'Bold');
        ylabel(['Recovered ', char(names(l))], 'FontWeight', 'Bold');
        filetype = '.fig';
        opt_3 = char(names(l));
        saveas(gcf, ['Local_PCA_n_clust_', opt_1, '_n_eig_', opt_2, '_var_', opt_3, filetype]);
    end

    % PLOT THE PCs SCORES FOR EACH LOCAL REGION
    
    for j = 1 : k
        figure
        subplot(1, 2, 1);
        plot(U_scores{j}(:, 1), U_scores{j}(:,2), 'b+');
        xlabel('1st U-score');
        ylabel('2nd U-score');
        subplot(1, 2, 2);
        plot(U_scores{j}(:, n_eig{j} - 1), U_scores{j}(:,n_eig{j}), 'b+');
        xlabel('last U-score');
        ylabel('2nd-last U-score');
        filetype = '.jpg';
        clust_number = num2str(j);
        opt_3 = char(clust_number);
        saveas(gcf, ['U_scores_n_clust_', opt_1, '_n_eig_', opt_2, '_clust_n.', opt_3, filetype]);
    end
end


fclose(fid);
%close all hidden
cd ..
                
    