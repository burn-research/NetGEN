function [y] = mole_to_mass(x, labels)
%  This function converts mole fractions to mass fractions of a given
%  species vector "labels", e.g. "labels = {'CH4', 'H2', ... }"


ns = length(labels);        % Number of species
np = length(x);             % Number of points

% Molecular weights
mw = mol_weight(labels);

for i = 1 : np

    M = 0;
    mi = zeros('like',mw);

    for j = 1 : ns
        M = M + x(i,j)*mw(j);
        mi(j) = x(i,j)*mw(j);
    end

    y(i,:) = mi/M;
end






end