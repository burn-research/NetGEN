function [y_axis, x_net, y_cfd, x_cfd] = plot_network_profile(idx, coord, var_name, var_cfd, path_to_output, opt)
% This function will plot the axial temperature profile of the CRN and
% compare it with the CFD solution. 
% INPUTS:
%   idx = clustering labelling vector
%   coord = matrix with the cartesian coordinates (mesh x-y-z)
%   T_cfd = vector of the temperature from CFD
%   path_to_output = path to the network output
%
% OUTPUTS:
%   y_axis = vector with the temperature profile from the network
%   x_net  = vector of axial coordinates to plot the network

% Find the cartesian coordinates that lie on the axis
id = [];            % Vector of indexes of points on the axis
x = coord(:,2);     % x-coordinate
y = coord(:,1);     % y-coordinate
clust_ax = [];      % Vector of clusters' indexes on the axis
for i = 1 : length(coord)

    if abs(x(i)) < 0.003
        id = [id; i];
        clust_ax = [clust_ax; idx(i)-1];
    end
end

y_axis = zeros('like', id);
if strcmp(var_name, 'T') == true
    for i = 1 : length(id)
        fname = append(path_to_output, 'Reactor.', num2str(clust_ax(i)), '/Output.out'); % Output name
        data = importdata(fname);
        val = data.data;
        y_axis(i) = val(5);
    end
else
    for i = 1 : length(id)
        fname = append(path_to_output, 'Reactor.', num2str(clust_ax(i)), '/Output.out'); % Output name
        data = importdata(fname);
        labels = data.textdata;
        val = data.data;
        for j = 1 : length(labels)
            ss = extractBefore(labels{j}, '_x');
            if strcmp(ss, var_name) == true
                idvar = j;
            end
        end
        y_axis(i) = val(idvar);
    end
end

[x_net, id_sort] = sort(y(id), 'ascend');
y_axis = y_axis(id_sort);
y_cfd = var_cfd(id(id_sort));
x_cfd = y(id(id_sort));

% Check if plot is true or false
if isfield(opt, 'Plot') == true 
    % Plot CFD and Network
    figure; hold on;
    [y_sort, id_sort] = sort(y(id), 'ascend');
    plot(y_sort, var_cfd(id(id_sort)),  'k-', 'LineWidth',1);
    plot(x_net, y_axis, 'r--', 'LineWidth',2);
    legend('CFD', 'CRN', 'interpreter', 'latex');
    ax = gca; ax.TickLabelInterpreter = 'latex';
end













end