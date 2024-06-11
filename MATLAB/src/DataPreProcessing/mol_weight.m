function [mw] = mol_weight(species)
% This function calculates the vector of the individual molecular weights
% of the given chemical species "species", where "species" is a vector of
% strings where the name of the species is contained, e.g. 
% species = {'CH4', 'H2', ... }

% Define atoms
Mc = 12;
Mh = 1;
Mn = 14;
Mo = 16;

% Initialize molecular weights
mw = zeros(length(species), 1);

for i = 1 : length(species)

    % Read the string
    sp = species{i};

    % Split the string
    sp_split = split(sp, '');
    sp_split = sp_split(2:end-1);

    % Initialize elements
    elem_list = {};     % Empty cell of elements strings
    elem_count = [];    % Array of elements stoich coefficients
    counter = 1;
    it = 1;
    conv = false;

    while conv == false
        
        % Iteration start
        if counter == 1
            elem_list{it} = sp_split{counter};
            if counter == length(sp_split)
                elem_count(it) = 1;
                conv = true;
            else
            
                % Check if next element in sp_split is a number or a letter
                if isempty(str2num(sp_split{counter+1})) == false % is a number
    
                    elem_count(it) = str2num(sp_split{counter+1});
                    counter = counter + 2;
    
                else % Is a letter, put 1 in the coeff
    
                    elem_count(it) = 1;
                    counter = counter + 1;
    
                end
    
                it = it + 1;
            end

        % If counter is not 1
        else 

            if counter > length(sp_split)
                conv = true;

            elseif counter == length(sp_split)
                elem_list{it} = sp_split{counter};
                elem_count(it) = 1;
                conv = true;

            else
                elem_list{it} = sp_split{counter};

                % Check if next element in sp_split is a number or a letter
                if isempty(str2num(sp_split{counter+1})) == false % is a number

                    elem_count(it) = str2num(sp_split{counter+1});
                    counter = counter + 2;
    
                else % Is a letter, put 1 in the coeff
    
                    elem_count(counter) = 1;
                    counter = counter + 1;

                end

                it = it + 1;
            end
        end
    end

    % Now get the molecular weight
    mi = 0;
    for j = 1 : length(elem_list)

        if strcmp(elem_list{j}, 'C') == true
            mi = mi + Mc*elem_count(j);
        elseif elem_list{j} == 'H'
            mi = mi + Mh*elem_count(j);
        elseif elem_list{j} == 'O'
            mi = mi + Mo*elem_count(j);
        elseif elem_list{j} == 'N'
            mi = mi + Mn*elem_count(j);
        end

    end

    mw(i) = mi;


end


















end