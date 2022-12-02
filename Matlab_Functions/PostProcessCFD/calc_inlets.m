function [inlets, bc_mass_flowrates] = calc_inlets(idx, inlet_streams, data_boundaries)
% This function, given a clustering label idx, will evaluate the correct
% inlets relative to inlet cluster, since a cluster can contain multiple
% inlets and/or a part of some inlet, so we need to create a function that
% will evaluate correctly the composition, the temperature and the mass
% flowrate of a stream entering the cluster
% 
% As input we have:
% idx = clustering labelling vector 
% inlet_streams = {stream_1, stream_2, ..., stream_n}
% data_boundaries = refer to function network_eval.m
% 
% As output we have:
% inlets = cell array with dimension nr x 1 where nr is the number of
%          reactors. The element in position j is a stream type array
%          containing:
%           inlets{j}.id = stream id (referring to Fluent zone id)
%           inlets{j}.T  = temperature [K]
%           inlets{j}.P = pressure [Pa]
%           inlets{j}.Y = cell type array {'H2:0.5', 'CH4:0.1', ...} mass
%                         composition
%           inlets{j}.Mf = mass flowrate [kg/s]
%
%

% Preliminary data
n_streams = length(inlet_streams);   % Number of streams in inlet
k = max(idx);                        % Number of cluster

stream_id = zeros(n_streams, 1);     % Vector with the id number of the streams
for j = 1 : n_streams
    stream_id(j) = inlet_streams{j}.id;
end

% Find the boundary cells
bc_cells = find(data_boundaries(:,3));   % Non-zero entry of the bc file

% Now create a matrix k x n_streams, where in the position (i,j) we have
% the mass flowrate of stream j entering in cluster i
% For the outlet does not matter, we compute only the mass flowrate
inlet_mass  = zeros(k, n_streams);
bc_mass_flowrates = zeros(k,2);

% Scan through the boundary cells to identify the inlets
for j = 1 : length(bc_cells)
    
    cell_id = data_boundaries(bc_cells(j), 1);  % Id of the cell
    zone_id = data_boundaries(bc_cells(j), 2);  % Id of the stream
    mi      = data_boundaries(bc_cells(j), 3);  % Mass flowrate in the cell
    
    clust_cell = idx(cell_id);                  % Cluster containing the inlet or outlet
    
    % If it's an outlet no problem
    if mi > 0
        bc_mass_flowrates(clust_cell,2) = bc_mass_flowrates(clust_cell,2) - mi;
    
    % If it's an inlet
    else
        ind = stream_id == zone_id;
        inlet_mass(clust_cell, ind) = inlet_mass(clust_cell, ind) - mi;
        bc_mass_flowrates(clust_cell,1) = bc_mass_flowrates(clust_cell,1) - mi;
        
    end
    
end

% Now we have the proportion of the streams entering each cluster, we have
% to mix them to find the right composition and temperature

% Final dictionary of the inlet streams associated with each cluster
inlets = cell(k, 1);

for j = 1 : k
    
    % If cluster k contains inlets
    if isempty(find(inlet_mass(j,:))) == false
        
        % For a certain cluster, if it is an inlet, we need to find which
        % streams are associated and in which proportion (mass flowrates of
        % each stream entering that cluster)
        ind_streams = find(inlet_mass(j,:) > 0);
        
        % If it's a multiple inlet
        if length(ind_streams) > 1
            disp('Multiple inlet detected');
            streams_to_mix = inlet_streams(ind_streams);  % Streams to be mixed
        
            % We need to calculate how much of each inlet streams go into
            % that cluster
            
            % Now scan through every inlet boundary cell belonging to the
            % cluster
            for i = 1 : length(streams_to_mix)
                
                % ID of the stream
                id_stream = streams_to_mix{i}.id;
                
                % Look for the id in the BC data
                ind = find(data_boundaries(:,2) == id_stream);
                
                % Now you have to take only the cells belonging to that
                % cluster
                cell_id = data_boundaries(ind, 1);
                mass_flow_stream = 0;
                
                % Scan through the cells belonging to the id of the stream,
                % locate the cells belonging to this cluster and calculate
                % the mass flowrate of the stream entering the cluster
                for l = 1 : length(cell_id)
                    clust_cell = idx(cell_id(l));
                    if clust_cell == j
                        mass_flow_stream = mass_flow_stream - data_boundaries(ind(l), 3);
                    end
                end 
                
                streams_to_mix{i}.Mf = mass_flow_stream;  % Modified mass flowrates
                
            end

            % Now we can mix the streams (refers to mix_streams.m)
            
            inlets{j} = mix_streams(streams_to_mix);
            
        % If it's a single inlet    
        elseif length(ind_streams) == 1
            
            inlets{j} = inlet_streams{ind_streams};
            
        end
        
    end
end
    
    
    


end

