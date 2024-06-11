function [output] = plot_clustering(coord, idx, dim)
% This function does a plot of the clustering in the geometrical space

set(0, 'defaulttextfontsize', 18);
set(0, 'defaultaxesfontsize', 18);
set(0, 'defaulttextinterpreter', 'latex');

k = max(idx); % Number of clusters

if dim == 2
    x = coord(:,1);
    y = coord(:,2);
    
    figure(1);
    scatter(y, x, 20, idx, 'filled'); colorbar; hold on;
    scatter(-y, x, 20, idx, 'filled');

    
    % Evaluate geometrical centroids
    C_geom = centroid_eval(coord, idx);
    
    % Display the number of the cluster of the scatter plot
    for i = 1 : length(C_geom)
        if C_geom(i,1) ~= 0 && C_geom(i,2) ~= 0
            text(C_geom(i,2), C_geom(i,1), num2str(i), 'Color', [1 0 0], 'FontSize', 8);
        end
    end
    
    xlabel('X coordinate [m]');
    ylabel('Y coordinate [m]');
    ax = gca; ax.TickLabelInterpreter = 'latex';

    cmap_name = append('parula(', num2str(k), ')');
    colormap(cmap_name);

    cb = colorbar;
    cb.Label.String = 'Cluster index';
    cb.Label.Interpreter = 'latex';


    elseif dim == 3
    x = coord(:,1);
    y = coord(:,2);
    z = coord(:,3);
    
    figure(1);
    scatter(x, z, 20, idx, 'filled'); colorbar;
    
    C_geom = centroid_eval(coord, idx);
    
    % Display the number of the cluster of the scatter plot
    for i = 1 : length(C_geom)
        if C_geom(i,1) ~= 0 && C_geom(i,3) ~= 0
            text(C_geom(i,1), C_geom(i,3), num2str(i), 'Color', [1 0 0], 'FontSize', 8);
        end
    end
end

output = true;