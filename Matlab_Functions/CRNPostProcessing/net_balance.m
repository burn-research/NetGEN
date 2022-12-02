function [global_imbalance, local_imbalance, mass_corrected] = net_balance(mass_flowrates, bc_mass_flowrates)
% This function evaluate the total mass imbalance of the network, both
% globally and locally (for each reactor)

% Preliminary information
[k, m] = size(mass_flowrates);
if k ~= m
    error('Connectivity matrix should be square');
end

% Initialize each reactor input and output
reactor_in  = zeros(k,1);
reactor_out = zeros(k,1);
local_imbalance = zeros(k,1);

% mass balance of the network
for j = 1 : k
    reactor_in(j)  = sum(mass_flowrates(:,j));          % Mass inflow from other reactors
    reactor_out(j) = sum(mass_flowrates(j,:));          % Mass outflow towards other reactors
    ext_input = bc_mass_flowrates(j,1);                 % Mass inflow from inlets
    ext_output = -bc_mass_flowrates(j,2);               % Mass inflow from outlets
    in = reactor_in(j) + ext_input;                     % Total input
    out = reactor_out(j) + ext_output;                  % Total output
    local_imbalance(j) = abs(in-out)/max(in,out);       % Local imbalance (relative %)
end

% Total input and total output in the network
total_in  = sum(reactor_in);
total_out = sum(reactor_out);

global_imbalance = abs(total_in - total_out)/max(total_in, total_out);

% Check if the mass balance is globally satisfied
mass_check = true;
if global_imbalance > 1e-12
    error('Check internal consistency of the network');
else
    
    % Check if the local mass balance is satisfied
    for j = 1 : k
        if local_imbalance(j) > 0.01 
            mass_check = false;
            str = append('Warning: mass imbalance over 1% detected in cluster ', num2str(j));
            disp(str);
        end
    end
end


if mass_check == false
    warning('Check local mass imbalance');
else
    disp('Mass balance is closed within the tolerance');
end

mass_corrected = mass_flowrates;

end

