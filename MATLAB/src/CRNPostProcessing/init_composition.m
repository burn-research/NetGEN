function [Y_list] = init_composition(data_all, idx, dim)
% This function initializes the composition of a reactor given the CFD
% solution and the idx labelling vector

% Check if labels are present
if isfield(data_all, 'Labels') == false
    error('Please specify labels CFD as field of the struct data_all');
end

% Exclude the first columns
if dim == 2
    start = 4;
elseif dim == 3
    start = 5;
end

comp = data_all.Solution(:,2:end);
labels_sp = data_all.Labels(start+1:end);

% Rewrite labels of species 
for i = 1 : length(labels_sp)
    li = split(labels_sp{i});
    if length(li) > 1
        LI = upper(li{2});
        labels_sp{i} = LI;
        % Check if contains molef
        if contains(LI, 'MOLE')
            LIS = extractAfter(LI, '-');
            labels_sp{i} = LIS;
        end
    end
end

% Cluster the data
comp_partition = clustering(comp, idx);

% Initialize the Y_list
k = max(idx);
Y_list = cell(k,1);
for i = 1 : k
    % Average composition
    comp_mean = mean(comp_partition{i});
    Yl = {};
    counter = 0;
    for j = 1 : size(comp_mean,2)
        if comp_mean(j) > 1e-6
            counter = counter + 1;
            ss = append(labels_sp{j}, ':', num2str(comp_mean(j)));
            Yl{counter} = ss;
        end
    end

    if isempty(Yl) == false
        Y_list{i} = Yl;
    else
        Y_list{i} = {'O2:0.21', 'N2:0.79'};
    end

end




end

