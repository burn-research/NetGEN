% LOCAL PCA - Performs local pca on a data set.
% 
% function [idx, uncentered_nz_X_k, F_k] = local_pca(X, F, names, C)
% INPUTS
% 
% X            = Matrix of state variables
% F            = Mixture fraction
% names        = State variables names (MUST be cellstr format)
% OPTIONAL:
% C            = cluster centroids initialization
%
% OUTPUTS
%
% idx               = Vector of indexes denoting the cluster number for
%                     each of the original data points
% uncentered_nz_X_k = Cell matrix of the partitioned data, i.e. each  cell
%                     is a cluster of points
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

function [idx, uncentered_nz_X_k, F_k, sq_rec_err, eigvec] = local_pca_lt_2(X, F, names, k, n_eigs, opt, C)

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
iter_max = 600;
eps_rec = 1.0;
eps_rec_min = 1.0e-02;
a_tol = 1.0e-16;
r_tol = 1.0e-08;

% Initial time
t0 = cputime;

% CENTER AND SCALE THE DATA
% Select centering and scaling criteria
% Before clustering we do not center, because  we remove the centroid
% inside the clustering we are centering (we do not center twice)

if isfield(opt, 'Center') == true
    cent = opt.Center;
    if cent ~= 0 || cent~= 1
        warning('Invalid centering criterion, data will be centerd by default');
        cent = 1;
    end
else
    warning('Centering option not specified, data will be centered by default');
end

pre_cent_crit = cent; %Refer to center.m for values
pre_scal_crit = 0; %Refer to scale.m for values

[cent_X, X_ave] = center(X, pre_cent_crit);
[scal_X, X_gamma] = scale(cent_X, X, pre_scal_crit);

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
    if nargin == 6
        
        % Option for initialization
        fprintf('\nClusters centroids not specified, initialization required \n');
        if isfield(opt, 'Initialization') == true
            init = opt.Initialization;
            if init < 1 || init > 3
                warning('Invalid initialization option selected, LPCA will be initialized uniform');
                init = 1;
            end
        else
            warning('Initialization option not selected explicitly. LPCA will be initialized uniformly');
            init = 1;
        end

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
        end
    else
        fprintf('\nClusters centroids initialized \n');
        opt_3 = 'Initialized';
    end
elseif VQ == 2
    opt_3 = 'Data partitioned on F';
    %The stoichiometric Z should be known.
end

% Initialization of eigenvectors
n_eigs_max = n_eigs;
eigvec = cell(k, 1);
gamma = cell(k, 1);

for j = 1 : k
    eigvec{j} = eye(columns, n_eigs);
    gamma{j} = ones(1, columns);
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

while ((convergence == 0) && (iter < iter_max))
    
    C_convergence = 0;
    eps_rec_convergence = 0;   
    fprintf('\nIteration n. %d, convergence %d \n', iter, convergence);  

    sq_rec_err = zeros(rows, k);
    
    if ((VQ == 2 || init == 3) && iter == 0) 
        F_min = min(F);
        F_max = max(F);
        F_stoich = input('\nInput the value of the stoichiometric mixture fraction: Fstoich = \n');
        %F_stoich = 0.4375;% DNS CO/H2 Flame
        [nz_X_k, nz_idx_clust] = condition(scal_X, F, k, F_min, F_max, F_stoich);
        C = zeros(k, columns);
        for j = 1 : k
            C(j, :) = mean(nz_X_k{j}, 1);
            [sort_eigval{j}, sort_eigvec{j}, eigval{j}, eigvec{j}, n_eig{j}, ...
                U_scores{j}, W_scores{j}, gamma{j}, scaled_data{j}, ...
                rec_data{j}] = pca_lt(nz_X_k{j}, cent_crit, scal_crit, 4, n_eigs);
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
    end
        
    [rec_err_min, idx] = min(sq_rec_err, [], 2);
    rec_err_min_rel = (rec_err_min);

    % Evaluate the global mean error
    eps_rec_new = mean(rec_err_min_rel);
    
    if VQ == 1
        % Partition the data into clusters
        [nz_X_k, nz_idx_clust, k] = partition(scal_X, idx, k);
    end
        
    fprintf('\nClusters dimension \n');
    disp(nz_X_k);
    
    % Evaluate the relative recontruction errors in each cluster
    rec_err_min_rel_k = cell(k, 1);
    for j = 1 : k
        rec_err_min_rel_k{j} = rec_err_min_rel(nz_idx_clust{j}, 1);
    end
    
    % Evaluate the mean error in each cluster
    eps_rec_new_clust = zeros(k, 1);
    size_clust = zeros(k, 1);
    for j = 1 : k
        eps_rec_new_clust(j) = mean(rec_err_min_rel_k{j});
        size_clust(j) = size(nz_X_k{j}, 1);
    end

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

        for j = 1 : k
            [sort_eigval{j}, sort_eigvec{j}, eigval{j}, eigvec{j}, n_eig{j}, ...
                U_scores{j}, W_scores{j}, gamma{j}, scaled_data{j}, rec_data{j}] ...
                = pca_lt(nz_X_k{j}, cent_crit, scal_crit, 4, n_eigs);
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

F_k = cell(k, 1);
unscaled_nz_X_k = cell(k, 1);
uncentered_nz_X_k = cell(k, 1);
for j = 1 : k
    F_k{j} = F(nz_idx_clust{j}, 1);
    unscaled_nz_X_k{j} = unscale(nz_X_k{j}, X_gamma);
    uncentered_nz_X_k{j} = uncenter(unscaled_nz_X_k{j}, X_ave(nz_idx_clust{j}, :));    
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
for j = 1 : k
    fprintf(fid, '\nCluster n. %d size = ', j);
    fprintf(fid, '%d \n', size_clust(j));
end
fprintf(fid, '\nFinal clusters centroids \n');
for j = 1 : k 
    for l = 1 : columns
    fprintf(fid,'%d \t', C(j, l));
    end
    fprintf(fid,'\n');
end

% SAVE DATA

save nz_X_k nz_X_k
save F_k F_k
save rec_X_data rec_X_data
save idx.out idx -ASCII -DOUBLE
save C.out C -ASCII -DOUBLE
save eigval eigval
save eigvec eigvec
save U_scores U_scores

% PLOTS

% PLOT OF THE ORIGINAL DATA VS RECONTRUCTED DATA

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
close all hidden
cd ..
                
    
    