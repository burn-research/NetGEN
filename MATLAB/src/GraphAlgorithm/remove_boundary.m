function [X_new, cell_to_remove, rows_to_remove] = remove_boundary(X)
rows_to_remove = [];
cell_to_remove = [];
count = 0;
[m,n] = size(X);
for i = 1 : m
    if X(i,2) == -1 || X(i,2) == 0
        count = count + 1;
        rows_to_remove(count,1) = i;
        cell_to_remove(count,1) = X(i,1);
    end
end

X_new = removerows(X, 'ind', rows_to_remove);

end

