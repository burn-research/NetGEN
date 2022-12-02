function [X_scaled] = log_scale(X)
% This function perform the logarithm of the data matrix X, avoiding
% negative logarithms

[m,n] = size(X);
X_scaled = zeros(m,n);

for i = 1 : m
    for j  = 1 : n
        if X(i,j) < 0
            X_scaled(i,j) = log(abs(X(i,j)));
        else
            X_scaled(i,j) = log(X(i,j));
        end
        
        if X_scaled(i,j) == -inf
            X_scaled(i,j) = -100.0;
        end
    end
end


            

end

