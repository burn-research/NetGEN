function [G] = weighted_graph(X, global_connect)
% 
% function [G] = weighted_graph(X, global_connect)
%
% X = data matrix
% global_coonect = m x 2 Matrix with the cell connections
% We want to build a weighted graph according to the euclidean distance
% between cells

[m,~] = size(global_connect);
weights = zeros(m,1);
for j = 1 : m
    weights(j) = sum((X(global_connect(j,1),:) - X(global_connect(j,2),:)).^2)/1e-6;
end


G = graph(global_connect(:,1)', global_connect(:,2)', weights', 'omitselfloops');

