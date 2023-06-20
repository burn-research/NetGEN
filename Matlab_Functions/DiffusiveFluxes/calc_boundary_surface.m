function [Ab] = calc_boundary_surface(G, subgraphs, idx, nodes, neighbor_data)
% This function will calculate the k x k surface boundary area between
% clusters given the global graph G, the subgraph dictionary subgraphs, the
% clustering labelling vector idx, the nodes vector nodes and the vector of
% areas

id1 = neighbor_data(:,4);
id2 = neighbor_data(:,5);
areas = neighbor_data(:,7);

% Check that subgraph number is the same as clusters
k = max(idx);
if k ~= length(subgraphs)
    error('Number of cluster different from length of subgraphs');
end

% For each subgraph, we need to find the neighbors of each point in the
% global graph G. Points that are in a different subgraph (clusters) will
% be counted and then the area will be calculated
Ab = zeros(k, k);
for i = 1 : k
    for j = 1 : k

        if i ~= j

            Hi = subgraphs{i};
            Hj = subgraphs{j};
    
            % Extract nodes
            Ni = nodes(idx==i);
            Nj = nodes(idx==j);
    
            % Neighbor nodes of Hi in subgraph Hj
            next = [];
            ai = [];
            count = 0;
    
            for l = 1 : length(Ni)
    
                % Find neighbors in G
                nn = neighbors(G, Ni(l));   % Neighbors of the nodes in Hi
    
                % Check if they are in cluster j
                cext = idx(nn); % External cluster but we have to take only j

                for s = 1 : length(cext)
                    if cext(s) == j
                        count = count + 1;
                        next(count,1) = nn(s);
                        next(count,2) = Ni(l);

                        % Find corresponding area
                        ida = find(id2 == Ni(l) & id1 == nn(s));
                        if isempty(ida)
                            ida = find(id1 == Ni(l) & id2 == nn(s));
                        end
                        if isempty(ida)
                            error('Cannot find neighbros');
                        end
                        ai(count) = areas(ida);
                    end
                end
            end
    
            Ab(i,j) = sum(ai);

        end

    end
end




end

