function [f] = mixture_fraction(X, labels, fuel)
% This function will calculate the mixture fraction according to the
% Bilger's definition. 
% INPUTS:
%   X = m x ns matrix with mass fractions of the species
%   labels = ns x 1 cell array with the name of the species
%   fuel = stream-like cell array (see network_eval.m) with fuel
%   oxidizer = stream-like cell array (see network_eval.m) with oxidizer
%
% OUTPUTS:
%   f = m x 1 vector with the pointwise mixture fraction

[m,ns] = size(X);
if ns ~= length(labels)
    error('Size of species matrix and labels not matching');
end

% Define parameters
MW_H = 1.0;
MW_O = 16.0;
MW_C = 12.0;
MW_N = 14.0;

% We need to identify for each species, and so for each column, how much
% carbon, oxygen, hydrogen or nitrogen there is
C_col = zeros(ns,1)';
O_col = zeros(ns, 1)';
H_col = zeros(ns, 1)';
N_col = zeros(ns, 1)';

for i = 1 : ns

    s = split(labels{i});
    sp_name = extractAfter(s{2}, '-');
    ch = split(sp_name, '');

    % Initialize element list as an empty cell
    elem_list = {};
    elem_count = [];
    count = 1;

    % This is a char vector with the name of the molecule, e.g. 
    % ch = 'c', '2', 'h', '3'
    for j = 2 : length(ch) - 1 % Exclude first and last because it is empty

        if j == 2 && count == 1
            elem_list{count} = ch{j};

            % Check if the next element of the vector is a
            % letter
            if isnan(str2double(ch{j+1})) == true
                elem_count(count) = 1;
            else
                elem_count(count) = str2double(ch{j+1});
            end

            count = count + 1;

        % In this case ch{j} can be a letter or a number, if it is a number
        % we skip cause it will be considered here
        elseif j > 2 && isnan(str2double(ch{j})) == true

            elem_list{count} = ch{j};

            % Check if the next element of the vector is a
            % letter
            if isnan(str2double(ch{j+1})) == true
                elem_count(count) = 1;
            else
                elem_count(count) = str2double(ch{j+1});
            end

            count = count + 1;

        end


    end

    % Check if elem_count and elem_list have the same length
    if length(elem_list) ~= length(elem_count)
        error('List of elements have different dimension of list of stoichiometric coefficients');
    end

    % Calculate the mass of each element present
    mc = 0;
    mn = 0;
    mo = 0;
    mh = 0;
    
    for j = 1 : length(elem_list)
        if strcmp(elem_list{j},'c') == true || strcmp(elem_list{j},'C') == true
            mc = elem_count(j)*12;
        elseif strcmp(elem_list{j},'n') == true || strcmp(elem_list{j},'N') == true
            mn = elem_count(j)*14;
        elseif strcmp(elem_list{j},'o') == true || strcmp(elem_list{j},'O') == true
            mo = elem_count(j)*16;
        elseif strcmp(elem_list{j}, 'h') == true || strcmp(elem_list{j},'H') == true
            mh = elem_count(j)*1;
        end
    end

    mtot = mc + mn + mo + mh;

    C_col(i) = mc/mtot;
    H_col(i) = mh/mtot;
    O_col(i) = mo/mtot;
    N_col(i) = mn/mtot;

end

% Now we need to calculate the mass fractions of C, O and H in the fuel
fuel_comp = fuel.Y;
ns_f = length(fuel_comp);       % Number of species in the fuel
yf = zeros(ns_f,1)';
C_col_f = zeros(ns_f,1)';       % Fraction of C for each species of the fuel
H_col_f = zeros(ns_f,1)';       % Fraction of H for each species of the fuel

for i = 1 : ns_f

    ss = fuel_comp{i};                          % String in the format i.e. C2H6:0.1
    sp = extractBefore(ss,':');                 % Species name (i.e. C2H6)
    yf(i) = str2double(extractAfter(ss,':'));   % Mass fraction of the species

    ch = split(sp, '');

    % Initialize element list as an empty cell
    elem_list = {};
    elem_count = [];
    count = 1;

    % This is a char vector with the name of the molecule, e.g. 
    % ch = 'c', '2', 'h', '3'
    for j = 2 : length(ch) - 1 % Exclude first and last because it is empty

        if j == 2 && count == 1
            elem_list{count} = ch{j};

            % Check if the next element of the vector is a
            % letter
            if isnan(str2double(ch{j+1})) == true
                elem_count(count) = 1;
            else
                elem_count(count) = str2double(ch{j+1});
            end

            count = count + 1;

        % In this case ch{j} can be a letter or a number, if it is a number
        % we skip cause it will be considered here
        elseif j > 2 && isnan(str2double(ch{j})) == true

            elem_list{count} = ch{j};

            % Check if the next element of the vector is a
            % letter
            if isnan(str2double(ch{j+1})) == true
                elem_count(count) = 1;
            else
                elem_count(count) = str2double(ch{j+1});
            end

            count = count + 1;

        end


    end

    % Calculate the mass of each element present
    mc = 0;
    mn = 0;
    mo = 0;
    mh = 0;
    
    for j = 1 : length(elem_list)
        if strcmp(elem_list{j},'c') == true || strcmp(elem_list{j},'C') == true
            mc = elem_count(j)*12;
        elseif strcmp(elem_list{j},'n') == true || strcmp(elem_list{j},'N') == true
            mn = elem_count(j)*14;
        elseif strcmp(elem_list{j},'o') == true || strcmp(elem_list{j},'O') == true
            mo = elem_count(j)*16;
        elseif strcmp(elem_list{j}, 'h') == true || strcmp(elem_list{j},'H') == true
            mh = elem_count(j)*1;
        end
    end

    mtot = mc + mn + mo + mh;

    C_col_f(i) = mc/mtot;
    H_col_f(i) = mh/mtot;

end

% Mass fraction of C and H in tne fuel
yc_f = sum(yf.*C_col_f);
yh_f = sum(yf.*H_col_f);

if yc_f + yh_f ~=1
    error('Sum of mass fractions in the fuel not zero');
end

% Now we need to calculate the mass fractions of O in the oxidizer
yo_o = 0.233;

% Now for each data point we can calculate the mixture fraction
f = zeros(m,1);
for i = 1 : m

    mc = sum(X(i,:).*C_col); 
    mo = sum(X(i,:).*O_col); 
    mh = sum(X(i,:).*H_col);
    mn = sum(X(i,:).*N_col);

    yc = mc/(mc + mo + mh + mn); 
    yo = mo/(mc + mo + mh + mn); 
    yh = mh/(mc + mo + mh + mn);

    %f(i) = (2*(yc)/12 + 2*(yh)/2 - (yo - yo_o)/16)/(2*yc_f/12 + yh_f/2 - yo_o/16);
    beta = 2*yc/12 + 0.5*yh - yo/16;
    betaf = 2*yc_f/12 + 0.5*yh_f;
    betao = -yo_o/16;

    f(i) = (beta-betao)/(betaf-betao);

end




