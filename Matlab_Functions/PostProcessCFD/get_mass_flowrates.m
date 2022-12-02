function [mass_flowrates] = get_mass_flowrates(alpha, bc_mass_flowrates)
% This function calculates the internal mass flowrates of a reactor network
% given the splitting ratio matrix and the external inputs and outputs

nr = size(alpha, 1);

A = eye(nr) - alpha';
b = bc_mass_flowrates(:,1);

M = linsolve(A,b);

mass_flowrates = zeros(nr, nr);
for i = 1 : length(M)
    mass_flowrates(i,:) = alpha(i,:) .* M(i);
end


end

