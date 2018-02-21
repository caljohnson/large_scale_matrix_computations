function [ y ] = row_stochastic_sparse_matrix_prod( n, Q, Jv, x )
%ROW_STOCHASTIC_SPARSE_MATRIX_PROD Matrix-vector product A^T x
%   where A is a sparse row-stochastic matrix:
%           A = Q + 1/n v e'
%       and x is a dense vector
%Inputs:  n  - dimension of matrix (A is n x n)
%         Q  - sparse directed graph matrix, before row-stochastic
%         Jv - vector of node indices j s.t. outdegree d_j = 0
%         x  - dense vector to mutliply by
%Outputs: y = A'x

y = Q'*x + ones(n,1)*sum(x(Jv))/n;

end

