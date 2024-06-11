function [inlet_id, inlet_mf] = get_mass_inlet(bc_inlet_data)
% This function returns the mass flowrates of the different inlets from the
% Fluent output file BC_mass_flowrates generated with the UDF

% Number of points in the data
n = length(bc_inlet_data);

% Initialize the counter
counter = 0;

% Initialize the inlets as empty list
inlet_id = [];
inlet_mf = [];

% Scan through the data
for i = 1 : n
    
    % If this is less than zero is an inlet
    if bc_inlet_data(i,3) < 0
        
        if counter == 0
            counter = counter + 1;
            inlet_id(counter) = bc_inlet_data(i,2);
            inlet_mf(counter) = -bc_inlet_data(i,3);
            
        else
            if isempty(find(inlet_id == bc_inlet_data(i,2)))
                counter = counter + 1;
                inlet_id(counter) = bc_inlet_data(i,2);
                inlet_mf(counter) = -bc_inlet_data(i,3);
                
            else
                ii = find(inlet_id == bc_inlet_data(i,2));
                inlet_mf(ii) = inlet_mf(ii) - bc_inlet_data(i,3);
            end
            
        end
    end
end



end

