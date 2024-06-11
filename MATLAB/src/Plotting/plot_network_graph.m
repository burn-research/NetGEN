function [output] = plot_network_graph(mass_flowrates, C, net_volumes, coord, idx)
% This function will plot a directed graph representing a certain reactor
% network

% Check network dimensions
[m,n] = size(mass_flowrates);
if m ~= n
    error('Mass flowrates matrix should be squared');
end

% Check the geometric dimension of the problem
[~, dim] = size(C);

% Log of the size of the reactors
r_size = log10(net_volumes).^2;
for i = 1 : length(r_size)
    if r_size(i) < 0
        r_size(i) = 0.1;
    end
end

mass_relative = 3*mass_flowrates./max(mass_flowrates);
h = digraph(mass_relative, 'omitselfloops');

% Create the digraph of the network
switch dim
    case 2
        figure;
        plot(h, 'XData', C(:,1), 'YData', C(:,2), 'MarkerSize', r_size, 'LineWidth', h.Edges);
    case 3
        figure;
        hold on;
        scatter(coord(:,1), coord(:,3), 10, idx, 'filled');
        
        plot(h, 'XData', C(:,1), 'YData', C(:,3), 'MarkerSize', r_size, 'NodeFontSize', r_size, 'NodeColor', 'r', 'EdgeColor', 'r', ...
            'NodeLabelColor', 'r', 'LineWidth', h.Edges.Weight);
        
end


output = h;


end

