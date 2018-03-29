function [ d_X_X_i ] = distance_matrix( X )
%Input: Data points X
%Output: Distance matrix d(x,x_i), where
%       d(i,j) = ||x_i - x_j||_2^2
N = size(X, 2);
d_X_X_i = zeros(N,N);
for i = 1:N
    for j = 1:N
        d_X_X_i(i,j) = mydistance(X(:,i),X(:,j));
    end
end

end

