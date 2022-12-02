function [x, y, yi, xi, data_q] = plot_contours(data, coord)
% This function will plot contours of the vector data on the 2D physical
% coordinates coord

x = coord(:,1); y = coord(:,2);

% Create the grid
xx = linspace(0, max(x), 1000);
yy = linspace(0, max(y), 1000);

[xi,yi] = meshgrid(xx,yy);

data_q = griddata(x, y, data, xi, yi);

plotting = 'no';
switch plotting
    case 'yes'
        plt = contourf(xi, yi, data_q, 200, 'LineStyle','none');
        xlabel('x [m]', 'Interpreter', 'latex'); ylabel('y [m]', 'Interpreter', 'latex');
        title('Contour', 'Interpreter', 'latex');
        
        ax = gca; ax.TickLabelInterpreter = 'latex';
        
        % Label the colorbar
        cb = colorbar;
        cb.TickLabelInterpreter = 'latex';
        
        colormap('jet(50)');
    case 'no'
        plt = true;
end


end

