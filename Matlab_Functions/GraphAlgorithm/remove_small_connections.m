function [mass_flow_new] = remove_small_connections(mass_flowrates)
% This function will remove small connections in reactor networks

[n,m] = size(mass_flowrates);
if n ~= m
    error('Dimension of the network not consistent');
end

% Initialize the new mass flowrates matrix
mass_flow_new = mass_flowrates;

% Check the connections in terms of relative amounts
for i = 1 : n
        Mi = sum(mass_flowrates(i,:));      % Outflow of reactor i
        alphai = mass_flowrates(i,:)/Mi;    % Relative outflow
        for j = 1 : n
            if alphai(j) < 1e-4
                mess = append('Small connection between reactor ', num2str(i), ' and ', num2str(j), ' it will be removed');
                disp(mess);
                mass_flow_new(i,j) = 0;
            end
        end
end
                
        


end

