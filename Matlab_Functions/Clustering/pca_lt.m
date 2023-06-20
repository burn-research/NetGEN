% Main function for PCA. The function provides a tool for data processing
% (centering and scaling), PCA analysis, PCs retention and original data
% reconstruction
%
% function [sort_eigval, sort_eigvec, ret_eigval, ret_eigvec, n_eig, U_scores, W_scores, gamma, scaled_data, rec_data, X_ave] = pca_lt(raw_data, cent_crit, scal_crit, stop_rule, inputs)
%


function [sort_eigval, sort_eigvec, ret_eigval, ret_eigvec, n_eig, U_scores, W_scores, gamma, scaled_data, rec_data, X_ave] = pca_lt(raw_data, cent_crit, scal_crit, stop_rule, inputs)

% CHECKING THE INPUTS

dim = size(raw_data);

if length(dim) ~= 2
    error('Matrix must be two-dimensional')
end

if nargin == 0
    error('You must provide the data matrix X \n');
end

if nargin == 1
    cent_crit = input('Choose the centering criterion: 0) No 1) Single 2) Double: \n');
    scal_crit = input('Scaling criterion: 0) No 1) STD 2) RANGE 3) PARETO 4) VAST 5) LEVEL 6) MAX  \n');
    stop_rule = input('Select the criterium to retain the PCs: 1)% Variance 2)PCs size 3)Broken Stick 4)Scree graph 5)Equality test \n');
    if stop_rule == 1
        var_lev = input('Select the variance to preserve = ');
    elseif stop_rule == 2
        threshold_rule = input('Select the threshold for PCs retention: 1) Kaiser rule 2) Joliffe rule \n');  
    elseif stop_rule == 4    
        n_eig = input('Choose the number of eigenvalues to retain ');  
    elseif stop_rule == 5
        alpha = input('Select the confidence level for the critical chi-square value: alpha = [0, 1] \n');
    end
end

if nargin == 2
    scal_crit = input('Scaling criterion: 0) No 1) STD 2) RANGE 3) PARETO 4) VAST 5) LEVEL 6) MAX  \n');
    stop_rule = input('Select the criterium to retain the PCs: 1)% Variance 2)PCs size 3)Broken Stick 4)Scree graph 5)Equality test \n');
    if stop_rule == 1
        var_lev = input('select the percentage of variance you want to preserve = ');
    elseif stop_rule == 2
        threshold_rule = input('Select the threshold for the PCs retention: 1) Kaiser rule 2) Joliffe rule \n');   
    elseif stop_rule == 4    
        n_eig = input('Choose the number of eigenvalues to retain ');   
    elseif stop_rule == 5
        alpha = input('Select the confidence level used to find the critical chi-square value for the significance test: alpha = [0, 1] \n');
    end
end

if nargin == 3
    stop_rule = input('Select the criterium to retain the PCs: 1)% Variance 2)PCs size 3)Broken Stick 4)Scree graph 5)Equality test \n');
    if stop_rule == 1
        var_lev = input('select the percentage of variance you want to preserve = ');
    elseif stop_rule == 2
        threshold_rule = input('Select the threshold for the PCs retention: 1) Kaiser rule 2) Joliffe rule \n');    
    elseif stop_rule == 4    
        n_eig = input('Choose the number of eigenvalues to retain ');   
    elseif stop_rule == 5
        alpha = input('Select the confidence level used to find the critical chi-square value for the significance test: alpha = [0, 1] \n');
    end
end

if nargin == 4
    if stop_rule == 1
        var_lev = input('select the percentage of variance you want to preserve = ');
    elseif stop_rule == 2
        threshold_rule = input('Select the threshold for the PCs retention: 1) Kaiser rule 2) Joliffe rule \n');    
    elseif stop_rule == 4    
        n_eig = input('Choose the number of eigenvalues to retain ');   
    elseif stop_rule == 5
        alpha = input('Select the confidence level used to find the critical chi-square value for the significance test: alpha = [0, 1] \n');
    end
end

if nargin == 5
    if stop_rule == 1
       var_lev = inputs;
    elseif stop_rule == 2
        threshold_rule = inputs; 
    elseif stop_rule == 4
        n_eig = inputs;
    elseif stop_rule == 5
        alpha = inputs;
    end
end

% Definition of parameters

a_tol = 1.0e-10;

[rows columns] = size(raw_data);

% DATA PREPROCESSING: CENTERING

% To get the principal componenets the data shas to be processed in order
% to center the column of X (variables with zero mean). Then 
% we will evaluate the eigenvalues of the matrix 1/(n-1) * X'X, which is
% the sample covariance matrix of X. Beside classical centering other 2
% possibilities are available:
% 1) The column of X are left uncentered
% 2) X is double centered, i.e. both rows and columns have zero mean

[centered_data, X_ave] = center(raw_data, cent_crit);

% DATA PREPROCESSING: SCALING

% At this point we can scale the data. Without scaling PCA can be
% meaningless since the tecnique is size dependent and it will simply lead
% pick a principal component equal to the larger variable, i.e.
% temperature. Different choices are available:  
% 0) No scaling 
% 1) Auto-scaling (STD), each variable is normalized by its standard
% deviation 
% 2) RANGE each variable is normalized by its range 
% 3) PARETO, each variable is scaled by the suare root of its standard
% deviation  
% 4) VAST, each variable is scaled by the standard deviation and
% coefficient of variation 
% 5) LEVEL, each variable is normalized by the mean of the data

[scaled_data, gamma] = scale(centered_data, raw_data, scal_crit);

% EVALUATION OF THE SAMPLE COVARIANCE MATRIX

% We perform compute the covariance matrix of the data scaled and centered.

%cov_data = 1/rows * scaled_data' * scaled_data;
cov_data = cov(scaled_data);

lambda = eig(cov_data);
[eigenvectors, eigenvalues] = eig(cov_data); %#ok<NASGU>

% Now we order the eigenvalues in descendent order, from the highest to the
% lowest. We use the sort function which give us also the permutation
% indexes of the eigenvalues vector

[sort_eigval, sort_index] = sort(lambda, 'descend');

% We have to sort the eigenvectors according to the eigenvalues. We will
% exploit the sort indexes derived in the previous step.

[rows_eigenvectors columns_eigenvectors] = size(eigenvectors);

sort_eigvec = zeros(rows_eigenvectors, columns_eigenvectors);
for i = 1 : columns_eigenvectors
    sort_eigvec(:,i) = eigenvectors(:, sort_index(i));
end

% We have to check that the iegenvectors we found actually form an
% orthonormal basis for the data representation. We need to check, then,
% that the inner product ai*aj = delta_ij where delta_ij is the Kronecker
% delta and equals 1 if i=j and 0 otherwise.

transp_sort_eigvec = sort_eigvec';
%for i = 1 : columns_eigenvectors;
 %   for j = i+1 : columns_eigenvectors;
  %      if ((transp_sort_eigvec(i, :) * sort_eigvec(:,j)) > a_tol);
   %         transp_sort_eigvec(i, :) * sort_eigvec(:,j)
    %        error('The eigenvectors do not form an orthonormal basis');
     %   end
   % end
%end
%display('The eigenvectors form an orthonormal basis for the data');

for i = 1 : columns_eigenvectors;
    if ((transp_sort_eigvec(i, :) * sort_eigvec(:,i)) == 1);
        continue
    end
end
%display('The eigenvectors form a normal basis')


% PRINCIPAL COMPONENTS RETENTION

% Now we have to select the principal components that have to be retained
% to describe the data. The following guidelines have been proposed: 
% 1) TOTAL VARIANCE: We can retain the eigenvalues needed to account for a
% specific percentage of the total variance (i.e. 80%). The required number
% of PCs is then the smallest value of m for which this chosen percentage
% is exceeded. 
% 2) INDIVIDUAL VARIANCE: We can retain the components whose eigenvalues
% are greater than the average of the eigenvalues (Kaiser, 1960) or than
% 0.7 times he average of the eigenvalues (Joliffe 1972). For a correlation
% matrix this average equals 1. 
% 3) BROKEN STICK MODEL: We can select the retained PCs according to the
% Broken Stick Model
% 4) SCREE GRAPH: We can use the scree graph, a plot of the 
% eigenvalues agaist their indexes, and look for a natural break between
% the large and small eigenvalues. 
% 5) SIGNIFICANCE TEST: We can perform a test of significance of the larger
% eigenvalues. 

% STOP RULE 1 - CUMULATIVE PERCENTAGE OF TOTAL VARIATION

% Perhaps the most obvious criterion for choosing the PCs to be retained is
% to select a (cumulative) percentage of total variation which one desires
% that the selected PCs contribute, say 80% or 90%. The required number of
% PCs is then the smallest value of m for which this chosen percentage is
% exceeded. 

if stop_rule == 1
    tot_var = sum(sort_eigval);
    n_eig = 1;
    sum_var = 0;
    i = 1;

    while ((sum_var < var_lev) && (n_eig < columns_eigenvectors))
        sum_var = sum_var + sort_eigval(i)/tot_var;
        i = i + 1;
        n_eig = i - 1;
    end
    ret_eigval = sort_eigval(1:n_eig);
    ret_eigvec = zeros(rows_eigenvectors, n_eig);
    for i = 1 : n_eig;
        ret_eigvec(:, i) = sort_eigvec(:,i);
    end      
end


% STOP RULE 2 - Size of Variances of Principal Components

% The previous rule is equally valid whether a covariance or a correlation
% matrix is used to compute the PCs. The rule described in this section is
% constructed specifically for use with correlation matrices, although it
% can be adapted for some types of covariance matrices. The idea behind the
% rule is that if all elements of x are independent, then the PCs are the
% same as the original variables and all have unit variances in the case of
% a correlation matrix. Thus any PC with variance less than 1 contains less
% information than one of the original variables and so is not worth
% retaining. The rule, in its simplest form, is sometimes called Kaiser�s
% rule (Kaiser, 1960) and retains only those PCs whose variances lk exceed
% 1. It can be argued that a cut-off at lk = 1 retains too few variables.
% Consider a variable which, in the population, is more-or-less independent
% of all other variables. In a sample, such a variable will have small
% coefficients in (p - 1) of the PCs but will dominate one of the PCs,
% whose variance will be close to 1 when using the correlation matrix. As
% the variable provides independent information from the other variables it
% would be unwise to delete it. However, deletion will occur if Kaiser�s
% rule is used, and if, due to sampling variation, lk < 1. It is therefore
% advisable to choose a cut-off lower than 1, to allow for sampling
% variation. Jolliffe (1972) suggested, based on simulation studies, that
% 0.7 is roughly the correct level.

if stop_rule == 2
    n_eig = 1;
    if threshold_rule == 1
        eig_max = mean(sort_eigval);
        i = 1;
        if sort_eigval(:) > eig_max
            ret_eigvec = sort_eigvec;
            ret_eigval = sort_eigval;
            n_eig = size(sort_eigval, 1);
        else
            while ((sort_eigval(i) > eig_max) && (n_eig < columns_eigenvectors))
                i = i + 1;
                n_eig = i - 1;
            end
            ret_eigval = sort_eigval(1:n_eig);
            ret_eigvec = zeros(rows_eigenvectors, n_eig);
            for i = 1 : n_eig;
                ret_eigvec(:, i) = sort_eigvec(:,i);
            end
        end
    end
    
    if threshold_rule == 2
        eig_max = 0.6 * mean(sort_eigval);
        i = 1;
        if sort_eigval(:) > eig_max
            ret_eigvec = sort_eigvec;
            ret_eigval = sort_eigval;
            n_eig = size(sort_eigval, 1);
        else
            while ((sort_eigval(i) > eig_max) && (n_eig < columns_eigenvectors))
                i = i + 1;
                n_eig = i - 1;
            end
            ret_eigval = sort_eigval(1:n_eig);
            ret_eigvec = zeros(rows_eigenvectors, n_eig);
            for i = 1 : n_eig;
                ret_eigvec(:, i) = sort_eigvec(:,i);
            end
        end
    end
end

% STOP RULE 3 - BROKEN STICK MODEL

% An alternative way of looking at the sizes of individual variances is to
% use the so-called broken stick model. If we have a stick of unit length,
% broken into p segments, then it can be shown that the expected
% length of the kth longest segment is:
% 
% lambda* = sum(1/j)/p j = k:p
% 
% One way of deciding whether the proportion of variance accounted for by
% the kth PC is large enough for that component to be retained is to
% compare the proportion with lambda*. Principal components for which the
% proportion exceeds lambda* are then retained, and all other PCs deleted.

if stop_rule == 3
    n_eig = 1;
    i = 1;
    stick_stop =1;
    while ((stick_stop == 1) && (n_eig < columns_eigenvectors))
        broken_stick = 0;
        for j = i : rows_eigenvectors
            broken_stick = broken_stick + 1/j;
        end
        stick_stop = (sort_eigval(i) > broken_stick);
        i = i + 1;
        n_eig = i - 1;
    end
    ret_eigval = sort_eigval(1:n_eig);
    ret_eigvec = zeros(rows_eigenvectors, n_eig);
    for i = 1 : n_eig
        ret_eigvec(:, i) = sort_eigvec(:,i);
    end
end
    
% STOP RULE 4 - SCREE PLOT

% One way to determine the number of eigenvalues to retain is to look at
% the scree graph,which is simply the plot of the eigenvalues sorted in
% descrnding order against their indexes. The number of eigenvalues to
% retain is based on the observation of the index k at which the slopes of 
% lines joining the plotted points are �steep� to the left of k, and 
% �not steep� to the right of it. This value of k, defining an �elbow� 
% in the graph, is then taken to be the number of components to be
% retained. The formalization of this tecnique would require to look at the
% difference between l(k) - l(k-1) and to check when it becomes fairly 
% constant. This tecnique requires some degree of judgment and it is used
% in a mainly graphical way.

if stop_rule == 4
    ret_eigval = sort_eigval(1:n_eig);
    ret_eigvec = zeros(rows_eigenvectors, n_eig);
    for i = 1 : n_eig
        ret_eigvec(:, i) = sort_eigvec(:,i);
    end
end

% STOP RULE 5 - TEST OF SIGNIFICANCE OF THE LARGER EIGENVALUES

% The number of PCs to retain can be also decided  with a test of
% significance of the larger components. we test the hyphothesis that the
% last k population eigenvalues are small and equal. The implication is
% that the first sample components capture all the essential dimensions,
% whereas the last components reflect noise. H0k = yp-k+1 = yp-k+2 = ... =
% yp where y are the eigenvalues of the covariance matrix. If H0 is true
% the last components will tend to form a straight line with small slope in
% the scree graph plot. To test Hok we use a likelihood approach based on
% the average of the last eigenvalues of the covariance matrix.
%
% ave(lp-k+1) = sum(lp-k+1/k)
% 
% Then we test the statistic
% 
% u = (n-(2p+11)/6)*(k ln(ave(lp-k+1))-sum(lnlp-k+1))
%
% which has an approximate chi square distribution. We reject H0 if
% u>=chi_square(alpha, degrees of freedom). The parameter alpha is the
% confidence level parameter used to find the critical chi-square value. To
% carry out the procedure we start by testing H02. If the test is accepted
% we then test H03 and so on until a test is rejected for some value of k.
% The procedure tends to retain more PCs than are really necessary. For
% correlation matrices, Jolliffe (1970) found that the rule often
% corresponds roughly to choosing a cut-off of about 0.1 to 0.2, which is
% much smaller than both Kaiser and Joliffe criteria. Moreover, it is easy
% to construct examples where the method gives silly answers. For instance,
% if there is one near-constant relationship among the elements of x, with
% a much smaller variance than any other PC, then the procedure rejects
% H0,p-2 and declares that all PCs need to be retained, regardless of how
% nearly equal are the next few eigenvalues. 

if stop_rule == 5
    samples = rows;
    n_eig = 1;
    i = 1;
    rejection = 0;
    while ((rejection == 0) && (n_eig < columns_eigenvectors))
        sum_last_eigenvalues = 0;
        sum_ln_last_eigenvalues = 0;
        for j = i : length(sort_eigval)           
            sum_last_eigenvalues = sum_last_eigenvalues + sort_eigval(j);
            sum_ln_last_eigenvalues = sum_ln_last_eigenvalues + log(sort_eigval(j));
        end
        k = length(sort_eigval) - i + 1;
        ave_last_eigenvalues = sum_last_eigenvalues/k;      % Average of the last k eigenvalues    
        u = (samples - (2*length(sort_eigval)+11)/6)*(k*log(ave_last_eigenvalues)-sum_ln_last_eigenvalues);   % Test statistic
        df = 0.5*(k-1)*(k+2);                               % Degrees of freedom
        critical_chi_square = chi2inv(alpha, df);           % Critical chi square value for the test
        rejection = (u < critical_chi_square);              % Logical value for test rejection
        i = i + 1;
        n_eig = i - 1;
    end
    ret_eigval = sort_eigval(1:n_eig);
    ret_eigvec = zeros(rows_eigenvectors, n_eig);
    for i = 1 : n_eig;
        ret_eigvec(:, i) = sort_eigvec(:,i);
    end
end

% Stop rule = 6: hard thresholding on eigenvalues
if stop_rule == 6
    m = dim(1);             % Number of rows
    n = dim(2);             % Number of columns
    
    % Calculate beta
    if m > n
    	beta = n/m;
    else
        beta = m/n;
    end
    
    % Expression for omega
    omega = 0.56*beta^3 - 0.95*beta^2 + 1.82*beta + 1.43;
    
    % Get the singular values from the eigenvalues
    calc_eigs = 1;
    switch calc_eigs
        case 1
            sigma = zeros('like', eigenvalues);
            for i = 1 : length(eigenvalues)
                sigma(i) = (eigenvalues(i)*(rows-1))^0.5;
            end
        case 2
            [U,D,V] = svd(scaled_data, 'econ');
            sigma = diag(D);
    end
        
    
    % Calculate the coefficient from the matlab script
    coeff = optimal_SVHT_coef(beta, 'false');
    
    if abs((coeff-omega)/coeff) > 0.1
        warning('Coefficient and omega are different');
    end
    
    % Follow the script from the paper
    y = sigma;
    y( y < coeff * median(y) ) = 0; 
    
    % Calculate non zeros entries of y
    n_eig = length(find(y));
    
    if n_eig == 0
        n_eig = 1;
    end
    
    % Retain the eigenvectors
    ret_eigval = sort_eigval(1:n_eig);
    ret_eigvec = zeros(rows_eigenvectors, n_eig);
    for i = 1 : n_eig
        ret_eigvec(:, i) = sort_eigvec(:,i);
    end
end

% Stop rule == 7
% Based on the assumption that useless PC are only responsible of noise, so
% that the U scores distribution will be close to a normal distribution
if stop_rule == 7
    U_vec = sort_eigvec;
    U_scores = scaled_data * U_vec;
    
    n = length(sort_eigvec);
    % We need to assess the distribution of the U_scores, and decide which
    % has a normal distribution and then can be discarded
    
    for i = 1 : n
        Utest = U_scores(:,i);
        Uc = center(Utest, 1);
        Us = scale(Uc, Utest, 1);
        
        % Test the normal distribution of Us
%         alpha = inputs;
%         P = DagosPtest(Us, alpha);
%         disp(P);
%         if P >= alpha
%             break
%         end

        % Calculate the skewness of the distribution
        sk = skewness(Us);
        ku = kurtosis(Us);
        if abs(sk) < 0.3 && abs(ku) > 3 && i ~= 1
            break
        end
        
    end
    n_eig = i;
    
    % Retain the eigenvectors
    ret_eigval = sort_eigval(1:n_eig);
    ret_eigvec = zeros(rows_eigenvectors, n_eig);
    for i = 1 : n_eig
        ret_eigvec(:, i) = sort_eigvec(:,i);
    end
end

% GET THE PRINCIPAL COMPONENT SCORES

% Three different kinds of Principal Component scores can be employed:
%
% 1) U-scores = obtained by using the U-vectors, i.e. the eigenvectors of
%    the covariance matrix S. The resulting U-scores are uncorrelated and
%    have variances equalk to the corresponding latent roots.
% 2) W-scores = The U vectors are scaled by the inverse of the latent roots
%    square root, i.e. V = L^-0.5 * U. The W-scores are still uncorrelated
%    and have variances equal unity.

U_vec = ret_eigvec;
W_vec = zeros(columns, n_eig);
for j = 1 : n_eig
    W_vec(:, j) = (ret_eigvec(:, j) / ret_eigval(j)^0.5);
end

U_scores = scaled_data * U_vec;
W_scores = scaled_data * W_vec;

% GETTING THE OLD DATA BACK
%
% Now we want to get the old data back. We compute first the data without
% the eman and than we get the original row data by adding the mean: These
% data are still scaled (bounded between 0 and 1) and so we have to rescale
% them to get the original row data back.

rec_adj_data = U_scores * ret_eigvec';

% Now we have to rescale the data back, according to the scaling criterium
% adopted 

[rec_cent_data] = unscale(rec_adj_data, gamma);

% In the case when all the eigenvectors are preserved, we want to verufy
% that the recovered scaled data match the original scaled data

if (n_eig == columns)
    if ((centered_data(:,:) - rec_cent_data(:,:)) < a_tol)
    else
        disp('Error: the recovered scaled data are NOT equal to the raw scaled data');
    end
end

% save recovered_scaled_data.out rec_scal_data -ASCII

% According to the centering criterion adopted, we have to un-center the
% recovered data.

rec_data = uncenter(rec_cent_data, X_ave);

% In the case when all the eigenvectors are preserved, we want to verify
% that the recovered data match the original data. This is another check
% beside the previous one on scaled data.

if (n_eig == columns)
    if ((raw_data(:,:) - rec_data(:,:)) < a_tol)
    else
        disp('Error: the recovered data are NOT equal to the raw data');
    end
end

