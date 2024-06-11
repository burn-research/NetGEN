function [res] = sys(M, split_ratio, bc)
% This function express the mass balance for each reactor of the network as
% a function of M, which are the total output flowrates from each reactor,
% excluded the boundary conditions

res = M - split_ratio'*M - bc;
end

