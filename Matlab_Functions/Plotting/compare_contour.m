function [Y_net, Y_net_contour] = compare_contour(Y_cfd, var_name, path_to_output, idx, coord)
% This function plot the contour of the CFD vs the contour of the CRN
% solution

% Change folder
output_folder = append(path_to_output, '/');
cd(output_folder);

% Enter each reactor output and save the variable of interest
end_cycle = false;
ii = 0;

Y_net = [];

while end_cycle == false
    rname = append('Reactor.', num2str(ii));
    
    if isfolder(rname) == true
        cd(rname);
        
        % Import reactor data
        data_reactor = importdata('Output.out');
        labels = data_reactor.textdata;
        val = data_reactor.data;
        
        % Find variable of interest
        for j = 1 : length(labels)
            ss = extractBefore(labels{j}, '(');
            if strcmp(ss, var_name) == true
                id = j;
            end
        end
        Y_net(ii+1) = val(id);
    else
        end_cycle = true;
    end
    ii = ii + 1;
    cd ../
end

% Compare the contours
figure;
scatter(0.001+coord(:,2), coord(:,1), 10, Y_cfd, 'filled'); cb = colorbar;
cb.Label.String = var_name;

Y_net_contour = zeros('like', Y_cfd);
for i = 1 : length(idx)
    Y_net_contour(i) = Y_net(idx(i));
end
colormap("winter(20)")

% Get colorbar limits
ll = cb.Limits;

hold on;
scatter(-0.001-coord(:,2), coord(:,1), 10, Y_net_contour, 'filled'); cb = colorbar;
if strcmp(var_name, 'T[K](5)') == true
    var_name = 'T [K]';
end

cb.Label.String = var_name;
cb.TickLabelInterpreter = 'latex';
cb.Label.Interpreter = 'latex';
caxis(ll);
xlabel('X coordinate [m]');
ylabel('Y coordinate [m]');
ax = gca; ax.TickLabelInterpreter = 'latex';
fig = gcf; fig.Units = 'centimeters';
fig.Position = [15 15 14 16];
colormap('winter(20)');

hold on;
plot(zeros(100,1), linspace(min(coord(:,1)), max(coord(:,1)), 100), 'k--', 'LineWidth',2);

% Get also the error and plot it
err = zeros('like', Y_cfd);
for i = 1 : length(Y_cfd)
    err(i) = abs((Y_cfd(i) - Y_net(idx(i))));
end

text(-0.25, 0.65, '(a)', 'FontSize',18);
text(0.20, 0.65, '(b)', 'FontSize',18);

ax.XTickLabel = [];


Y_net_contour = Y_net_contour';

end

