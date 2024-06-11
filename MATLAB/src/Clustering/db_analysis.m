function [db, k_opt] = db_analysis(X, k_span, labels)
%
% function [db, k_opt] = db_analysis(X, k_span, labels)
%
% This function will perform an a-posteriori clustering evaluation using
% the metric DB index.

% Preliminary information
l = length(k_span);
db = zeros(l,1);

% Centering
cent_crit = input('Enter 1 for centering, 0 for no centering');


% Select scaling criteria
scal = input('Scaling criteria, available choises are: no, auto, range, pareto, vast, level, max, : ', 's');
if strcmp(scal, 'auto') == true
    scal_crit = 1;
elseif strcmp(scal, 'range') == true
    scal_crit = 2;
elseif strcmp(scal, 'pareto') == true
    scal_crit = 3;
elseif strcmp(scal, 'vast') == true
    scal_crit = 4;
elseif strcmp(scal, 'level') == true
    scal_crit = 5;
elseif strcmp(scal, 'max') == true
    scal_crit = 6;
elseif strcmp(scal, 'no') == true
    scal_crit = 0;
else
    disp('No criteria has been selected, auto scale by default');
    scal_crit = 1;
end

% Select the algorithm for clustering
alg = input('Enter the clustering algorithm, available choices are: k-means, gmm, lpca: ', 's');

% Data matrix, stopping criteria
switch alg
    case 'k-means'
        X_center = center(X, 1);
        X_scaled = scale(X_center, X, 1);
        
    case 'lpca'
        stop_crit = input('Stopping criteria for lpca: 1 for tot_var, 4 for fixed eigs: ');
        if stop_crit == 1
            input_lpca = input('Fraction of total variance to retain (from 0 to 1): ');
        elseif stop_crit == 4
            input_lpca = input('Number of eigenvectors to retain: ');
            if input_lpca > size(X, 2)
                error('Number of eigenvectors is higher than number of dimensions');
            end
        else
            fprintf('\n Stopping rule ot available, automatically set to 0.9 of total variance \n \n');
            stop_crit = 1;
            input_lpca = 0.9;
        end
        
    case 'gmm'
        X_center = center(X, 1);
    	X_scaled = scale(X_center, X, 1);
        options = statset('MaxIter',1000);
end

% Clustering evaluation
% Execute the clustering for each number of clusters in k_span
for j = 1 : l
    
    switch alg
        case 'k-means'
            idx = kmeans(X_scaled, k_span(j), 'Start', 'uniform', 'MaxIter', 1000);
            ev1 = evalclusters(X, idx, 'DaviesBouldin');
            db(j) = ev1.CriterionValues;
            
            
        case 'lpca'
   
            [idx, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, sq_rec_err] = ...
                local_pca_lt(X, labels, k_span(j), size(X,2), [], cent_crit, scal_crit, stop_crit, input_lpca);
            
            db(j) = db_pca(X, idx, sq_rec_err);
            
        case 'gmm'
            GMM = fitgmdist(X_scaled, k_span(j), 'Options', options, 'CovarianceType', 'diagonal', 'RegularizationValue', 0.05);
            idx = GMM.cluster(X_scaled);
            ev1 = evalclusters(X, idx, 'DaviesBouldin');
            db(j) = ev1.CriterionValues;
    end
    
end

% Plot results
figure;
plot(k_span, db, 'b--o', 'LineWidth', 3);
xlabel('Number of clusters');
ylabel('DB index');
tit = append('DB evaluation with ', alg);
title(tit);
ax = gca; ax.TickLabelInterpreter = 'latex';

[~, id_best] = min(db);
k_opt = k_span(id_best);

