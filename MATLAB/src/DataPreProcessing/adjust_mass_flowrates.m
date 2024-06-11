function [mass_adjusted] = adjust_mass_flowrates(mass_flowrates, bc_mass_flowrates)
%
% function [m_adjusted] = adjust_mass_flowrates(mass_flowrates, bc_mass_flowrates)
%

% Get mass splitting ratios
alpha = mass_split(mass_flowrates, bc_mass_flowrates);
for i = 1 : k
    for j = 1 : k
        if alpha(i,j) < 1e-5
            alpha(i,j) = 0;
        end
    end
end
        

% Solve the system to get the corrected mass flowrates
A  = eye(k) - alpha';
b  = bc_mass_flowrates(:,1);
mf = linsolve(A,b);
mass_adjusted = alpha.*mf;     

end

