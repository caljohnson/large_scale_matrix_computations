%HW 1 Problem 4d  - MAT226B WQ1028
%Carter Johnson

%use power method to compute dominant eigenvalue and eigenvector for
% the www0 matrix - using altered matrix-prod method -> A^T_alpha x

load('www0.mat')
load('x0.mat')

%create function to pass to power method
matrix_prod = @(x) row_stochastic_sparse_altered_matrix_prod(0.85, n, Q, Jv, x);

%use power method with parameters
eps = 10^(-12);
k_max = 10000;

%for initial condition x = e
e = ones(n,1);
[lambda1, x1] = power_method(matrix_prod, e, eps, k_max);

%for other initial condition x = x0 from 'x0.mat'
[lambda2, x2] = power_method(matrix_prod, x0, eps, k_max);

%compute relative eigenvalue residual
res1 = max(abs(matrix_prod(x1) - lambda1*x1))/max(abs(x1));
res2 = max(abs(matrix_prod(x2) - lambda2*x2))/max(abs(x2));

%compute 10 most important websites from PageRank
[PageRank1, I1] = sort(x1, 'descend');
[PageRank2, I2] = sort(x2, 'descend');
I1(1:10)
I2(1:10)