function [output] = post_process_results(path_to_output, path_to_input, opt)
% This function will post process the results of the CRN simulation by
% printing and plotting selected quantities
% INPUTS:
%   path_to_output = path to output folder 
%   path_to_input  = path to input.dic file in the form /path/to/intput/input.dic
%   opt = struct array with options
%
% OUTPUTS:
%   output = output flag

% Get outlet reactor
r_out = get_outlet_reactor(path_to_input);
if length(r_out) ~=1 
    r_out = r_out(1);
end

% Select species to export
if isfield(opt, 'ExportSpecies') == true
    sp = opt.ExportSpecies;
else
    sp = {'NO', 'NO2'};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% SPECIES EXTRACTION %%%%%%%%%%%%%%%%%%%%%%%%%%%
output.Species = zeros(length(sp),1);
output_file = append(path_to_output, '/Reactor.', num2str(r_out), '/Output.out');

% Import data from output file
data_output = importdata(output_file);
labels = data_output.textdata;
val = data_output.data;

% Search for the species to export
id_sp = [];
sp_strings = cell(length(sp),1);
for i = 1 : length(sp)
    sp_strings{i} = append(sp{i}, '_x');
end

% Check indexes of matching strings
count = 1;
for i = 1 : length(labels)

    % Extract string from labels
    ss = extractBefore(labels{i}, '(');
    
    % Check if strings compare
    for j = 1 : length(sp_strings)
        if strcmpi(sp_strings{j}, ss)
            id_sp(count) = i;
            count = count + 1;
        end
    end
end

% Check if all the species were found
if length(id_sp) < length(sp)
    error('Some species were not found in the output file');
elseif length(id_sp) > length(sp)
    error('Some redundant species was found');
end

% Extract the species from the output file
for i = 1 : length(sp)
    output.Species(i) = val(id_sp(i));
    output.SpeciesName{i} = sp{i};
end

end