function [mixed_stream] = mix_streams(streams)
% This function implement a mixer between multiple streams, calculating the
% new averaged composition and temperature. The input is a struct-type array
% containing multiple streams' information. In particular, every
% stream has the following attributes:
%
% stream.id = Identification number (refers to inlets in Fluent and networkEval.m)
% stream.T  = Temperature
% stream.P  = Pressure
% stream.Y  = Mass fractions in the form {'H2:0.5', 'CH4:0.5'}
% stream.Mf = Mass flowrate in kg/s
%
% The input is then a dictionary containing all the streams that has to be
% mixed
%
% The function will create a new stream with an averaged composition and
% temperature, in the same stream format

% Preliminary information
n_streams = length(streams);

% Initialize new stream
mixed_stream = struct;

% Calculate the new mass flowrate and new averaged temperature summing all the mass flowrates
% and evaluating the averaged temperature with mass flowrates

mdot_mix = 0;                       % Total mass flowrate
T_weight = 0;                       % sum(Ti*mdot_i) to be divided for mdot_mix
id_list = zeros(n_streams, 1);      % Vector of id's of the streams, to know which streams have been mixed

for j = 1 : n_streams
    id_list(j) = streams{j}.id;     % Update id list
    mdot = streams{j}.Mf;           % Mass flowrate of stream j
    T = streams{j}.T;               % Temperature of stream j
    mdot_mix = mdot_mix + mdot;     % Update the total mass flowrate
    T_weight = T_weight + T*mdot;   % Update temperature weighted sum
end

%% Averaged temperature
T_ave = T_weight/mdot_mix;

%% Evaluation of the averaged composition

% First of all, we create a global array with all the mass fractions
species_list = {};
masses = [];

% This counter will be updated every time a new species is found. If a
% species that already exist is found, then the mass flowrate of that
% species, intended mdot*yi is summed up in the vector species
count = 0;
for j = 1 : n_streams
    comp_i = streams{j}.Y;     % Composition vector of stream j
    mdot_i = streams{j}.Mf;    % Mass flowrate of stream j
    
    % Scan through the species in the composition vector comp i
    for i = 1 : length(comp_i)
        sp = extractBefore(comp_i{i}, ':');             % Species string
        yi = str2double(extractAfter(comp_i{i}, ':'));  % Mass fraction
        mi = mdot_i * yi;                               % Mass flowrate of species i
        
        % Scan through species_list if the species already exist
        % If species list is empty, just add the first species
        if count == 0
            species_list{1} = sp;
            masses(1,1) = mi;
            count = count + 1;
        
        % If it is not empty, then scan through already exisiting species
        else
            sp_exist = false;
            for l = 1 : count
                
                % If species are matched, add the masses
                if strcmp(sp, species_list{l}) == true
                    
                    sp_exist = true;                    % Update boolean variable
                    masses(l,1) = masses(l,1) + mi;     % The new mass flowrate of the i-th species
                    
                end
            end
            
            % If the species does not exist yet, add it to the list,
            % increasing the counter
            if sp_exist == false
                count = count + 1;
                species_list{count} = sp;
                masses(count,1) = mi;
            end
            
        end 
    end
end

% Now we can evaluate the new mass fractions by dividing the masses for the
% toal mass flowrate
yi_mix = masses/mdot_mix;
        
% We can round yi_mix for 5 significative numbers
yi_mix = round(yi_mix, 7);

% We can then make sure that the sum of the mass fraction is one by
% updating the last element of yi_mix as 1 - sum(y_mix(:,end-1))
yi_mix(end) = 1 - sum(yi_mix(1:end-1));

% Now we can rewrite the new composition dictionary in the stream format
n_species = count;
comp_mix = cell(count, 1);

for i = 1 : n_species
    comp_mix{i} = append(species_list{i}, ':', num2str(yi_mix(i)));
end

%% Evaluation of the pressure
% Take the minimum of the pressures
P0 = streams{1}.P;
for j = 2 : n_streams
    P = min(P0, streams{j}.P);
end

%% Now we can build the final mixed stream
mixed_stream.id = id_list;
mixed_stream.Y = comp_mix;
mixed_stream.T = T_ave;
mixed_stream.Mf = mdot_mix;
mixed_stream.P = P;


end

