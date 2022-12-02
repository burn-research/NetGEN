% This function provides a tool to center the data
% 
% INPUTS:
% 
% uncentered_data   = The matrix of uncentered variables
% cent_crit         = centering criterion
%
% OUTPUT
%
% centered_data     = Matrix of centered data
% X_ave             = Matrix of means to be subtracted

function [centered_data, X_ave] = center(uncentered_data, cent_crit)

% To get the principal componenets the data shas to be processed in order
% to center the column of X (variables with zero mean). Then 
% we will evaluate the eigenvalues of the matrix 1/(n-1) * X'X, which is
% the sample covariance matrix of X. Beside classical centering other 2
% possibilities are available:
% 1) The column of X are left uncentered
% 2) X is double centered, i.e. both rows and columns have zero mean

% Definition of parameters

[rows, columns] = size(uncentered_data);

mean_var = mean(uncentered_data, 1);
mean_obs = mean(uncentered_data, 2);
grand_mean = mean(mean_var);

% Relative tolerance

a_tol = 1e-06;

switch cent_crit
    case 0
        X_ave = zeros(rows, columns);
    case 1
        X_ave = repmat(mean_var, rows, 1);
    case 2
        X_ave = repmat(mean_var, rows, 1) + repmat(mean_obs, 1, columns) - repmat(grand_mean, rows, columns);
    otherwise
        error('Unknown centering criterion');
end

centered_data = uncentered_data - X_ave;

% We check now that each variable has zero mean

switch cent_crit
    case 1
        mean_vars_data = mean(centered_data, 1);
        for j = 1 : columns
            if (abs(mean_vars_data(j)) > a_tol)
                disp(mean_vars_data);
                error('The mean of each variable does not equal zero');
            end
        end
    case 2
        mean_vars_data = mean(centered_data, 1);
        mean_obs_data = mean(centered_data, 2);
        for j = 1 : columns
            if (abs(mean_vars_data(j)) > a_tol)
                error('The mean of each variable does not equal zero');
            end
        end
        for i = 1 : rows
            if (abs(mean_obs_data(i)) > a_tol)
                error('The mean of each row does not equal zero');
            end
        end
end

%save centered_data centered_data
%save X_ave X_ave