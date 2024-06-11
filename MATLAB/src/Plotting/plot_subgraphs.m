function [H] = plot_subgraphs(cell_id, idx, G, coord)
% This function will plot the graphs associated with each cluster, in order
% to locate disconnected graph for a detailed analysis

% Partition the cells  and the coordinates according to the labels vector idx
cell_clust = clustering(cell_id, idx);
coord_clust = clustering(coord, idx);

% Number of clusters
k = max(idx);

% Check the geometrical dimension (2D ord 3D)
dim = size(coord, 2);

% Scan through each cluster
switch dim
    case 2
        for j = 1 : k
            cell_clust_j = cell_clust{j};
            coord_clust_j = coord_clust{j};
            H = subgraph(G, cell_clust_j);
            bins = conncomp(H, 'OutputForm', 'cell');
            figure;
            p = plot(H, 'XData', coord_clust_j(:,1), 'YData', coord_clust_j(:,2));
            tit = append('Cluster ', num2str(j));
            title(tit);
            
            % Sort the subgraphs for size
            binsize = cellfun('length', bins);
            [binsize, ind] = sort(binsize, 'descend');
            bins = bins(ind);
            
            for i = 2 : length(bins)
                highlight(p, bins{i}, 'NodeColor', 'r', 'Marker', 'x');
            end
        end
        
    case 3
        for j = 1 : k
            cell_clust_j = cell_clust{j};
            coord_clust_j = coord_clust{j};
            H = subgraph(G, cell_clust_j);
            bins = conncomp(H, 'OutputForm', 'cell');
            figure;
            p = plot(H, 'XData', coord_clust_j(:,1), 'YData', coord_clust_j(:,2), 'ZData', coord_clust_j(:,3));
            tit = append('Cluster ', num2str(j));
            title(tit);
            % Sort the subgraphs for size
            binsize = cellfun('length', bins);
            [binsize, ind] = sort(binsize, 'descend');
            bins = bins(ind);
            
            for i = 2 : length(bins)
                highlight(p, bins{i}, 'NodeColor', 'r', 'Marker', 'x');
            end
        end
end
        
    
    


end

