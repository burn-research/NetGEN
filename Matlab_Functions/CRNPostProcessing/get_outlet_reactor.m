function [out] = get_outlet_reactor(filename)
% This function read the input.dic dictionary of the reactor network and
% get the outlet reactor. Filename should be passed as a string

lines = readlines(filename);

% Scan through the lines and detect the one starting with the desired
% keyword
for i = 1 : length(lines)
    l = split(lines(i));
    if strcmp(l(1), '@OutputStreams') 
        out = str2num(l(2));
    end
end




end

