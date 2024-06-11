function [Y_net] = plot_network_contour(idx, coord, var_name, path_to_output)
% This function plot the contour obtained from the network
% INPUTS
%   idx = clusteruing label vector
%   coord = coordinate of cell centroids from CFD simulations
%   var name = string corresponding to T or species name (e.g. NO)
%   path_to_output = path pointing to Output/ folder of the desired case
%
% OUTPUTS
%   output = true if everything was alright

% Initialize output vector
Y_net = zeros(length(idx),1);
if strcmp(var_name, 'idx') == true
    Y_net = idx;
    figure;
    pp = plot_contours(Y_net, coord);
    cm = append('parula(', num2str(max(idx)), ')');
    colormap(cm);
    hold on;
    id_inj = coord(:,1) < 0;
    scatter(coord(id_inj,2), coord(id_inj,1), 20, idx(id_inj), 'filled');
else

% Create data partitions
idx_clust = clustering([1:1:length(Y_net)]', idx);

% Initialize column index corresponding to var name
ind = 0;
for i = 1 : length(idx_clust)

    ri = i-1;
    rname = append(path_to_output, 'Reactor.', num2str(ri), '/Output.out');
    data = importdata(rname);
    
    % At the beginning we need to find the column corresponding to
    % var_name
    if i == 1
        labels = data.textdata;
        for j = 1 : length(labels)

            s = labels{j};
            si = extractBefore(s, '_x');

            if strcmp(si, var_name) == true
                ind = j;
            end
        end

        if strcmp(var_name, 'T') == true
            ind = 5;
        end
    end

    val = data.data;
    Y_net(idx_clust{i}) = val(ind);
end

% Plot results
plt = 'no';
switch plt
    case 'scatter'
        figure;
        scatter(coord(:,2), coord(:,1), 10, Y_net, 'filled');
        hold on;
        scatter(-coord(:,2), coord(:,1), 10, Y_net, 'filled');
        colormap('winter(20)');
        cb = colorbar;
        cb.Label.String = var_name;
        ax = gca; ax.TickLabelInterpreter = 'latex';
    case 'contour'
        figure;
        pp = plot_contours(Y_net, coord);
    case 'mix'
        figure;
        pp = plot_contours(Y_net, coord);
        hold on;
        id_inj = coord(:,1) < 0;
        scatter(coord(id_inj,2), coord(id_inj,1), 20, idx(id_inj), 'filled');
end
end

        








end

