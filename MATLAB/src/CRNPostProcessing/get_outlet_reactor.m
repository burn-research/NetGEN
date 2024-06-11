function [out] = get_outlet_reactor(filename)
% This function read the input.dic dictionary of the reactor network and
% get the outlet reactor. Filename should be passed as a string

lines = readlines(filename);

% Scan through the lines and detect the one starting with the desired
% keyword
out = [];
for i = 1 : length(lines)
    l = split(lines(i));
    if strcmp(l(1), '@OutputStreams') 
        out(1) = str2num(l(2));

        % Check the next lines, if they are a number they are also outputs
        stop = false;
        counter = 1;
        while stop == false
            ll = split(lines(i + counter));
            a = str2num(ll(1));
            if isnumeric(a) == true
                out(counter+1) = a;
            else
                stop = true;
                break
            end
            counter = counter + 1;
        end
    end
end




end

