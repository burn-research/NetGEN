function [mass_flow_corrected, mass_imbalance] = mass_balance(mass_flowrates, bc_flows)
% This function will compute the new mass flowrates to satisfy the mass
% balance of the network

% Number of reactors
k = size(mass_flowrates, 1);

%% Remove very small mass flowrates
for i = 1 : k
    for j = 1 : k
        if mass_flowrates(i,j) < 1e-5
            mass_flowrates(i,j) = 0;
        else
            mass_flowrates(i,j) = round(mass_flowrates(i,j), 6, 'significant');
        end
    end
end

%% Calculate split ratios
% Remember that mass flowrates (i,j) is the mass flowrate from reactor i to
% reactor j

% Calculate the mass flow splits
% Flow_split(i,j) = fraction of mass from reactor i to j
flow_split = zeros(k, k);
for i = 1 : k
    flow_split(i, :) = mass_flowrates(i,:)/sum(mass_flowrates(i,:));
end

alpha = flow_split;

%% Calculate total outflows
M = sum(mass_flowrates, 2);   % Vector of the outflows

%% Define the system
    function [res] = reac_balance(dM)
        res = zeros(k,1);
        for l = 1 : k
            res(l) = (M(l) + dM(l)) - alpha(:,l)'*(M + dM) - bc_flows(l,1) - bc_flows(l,2);
        end
    end

%% Linear system
iter = 1;
convergence = false;
itmax = 100;
A = (eye(k)-alpha');
while convergence == false && iter < itmax
    if iter ~= 1
        M = M + x;
    end
    b = bc_flows(:,1) + bc_flows(:,2) - (eye(k) - alpha')*M;
    x = lsqminnorm(A,b);
    resid = norm(A*x - b);
    if resid == 0
        convergence = true;
        disp('Convergence riched, mass balance perfectly matching');
    else
        iter = iter + 1;
    end
    
    if iter == itmax
        disp('Residuals decreased');
    end
end

%% New values of mass flowrates
mass_flow_corrected = zeros(k,k);
for i = 1 : k
    for j = 1 : k
        mass_flow_corrected(i,j) = alpha(i,j)*M(i);
    end
end

mass_imbalance = reac_balance(x);

end

