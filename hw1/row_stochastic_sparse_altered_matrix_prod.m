function [ y ] = row_stochastic_sparse_altered_matrix_prod( alpha, n, Q, Jv, x )
%ROW_STOCHASTIC_SPARSE_ALTERED_MATRIX_PROD Matrix-vector product A_alpha^T x
%   where A is a sparse row-stochastic matrix:
%           A = Q + 1/n v e'
%   and A_alpha is altered:
%           A_alpha = alpha A + (1-alpha)e e^T /n            
%       and x is a dense vector
%Inputs:  alpha - alteration factor, 0 <= alpha <= 1
%         n  - dimension of matrix (A is n x n)
%         Q  - sparse directed graph matrix, before row-stochastic
%         Jv - vector of node indices j s.t. outdegree d_j = 0
%         x  - dense vector to mutliply by
%Outputs: y = A_alpha^T x

y = alpha*Q'*x;
y = y + (alpha/n)*ones(n,1)*sum(x(Jv)) + ((1-alpha)/n)*ones(n,1)*sum(x);
end

