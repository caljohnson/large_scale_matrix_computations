%HW 1 Problem 4c  - MAT226B WQ1028
%Carter Johnson

%use power method to compute dominant eigenvalue and eigenvector for
% the matrix from problem 2

%matrix from problem 2
A = [ 0 0 0 0 1/2 0 1/2 0 0 0; 
0 1/4 0 0 0 1/4 1/4 0 1/4 0; 
1/10 1/10 1/10 1/10 1/10 1/10 1/10 1/10 1/10 1/10;
0 0 0 1/2 0 0 1/2 0 0 0;
0 0 0 0 0 0 1/3 0 1/3 1/3;
1/10 1/10 1/10 1/10 1/10 1/10 1/10 1/10 1/10 1/10;
1/10 1/10 1/10 1/10 1/10 1/10 1/10 1/10 1/10 1/10;
0 0 1/5 1/5 0 0 1/5 1/5 1/5 0;
0 1/5 0 1/5 0 0 1/5 1/5 0 1/5;
1/6 0 1/6 1/6 0 0 1/6 0 1/6 1/6;];

%get sparse matrix Q, before we added the 1/10 rows
n=10;
A(A<=1/10)=0;
Q = sparse(A);
Jv = [3 6 7];

%create function to pass to power method
matrix_prod = @(x) row_stochastic_sparse_matrix_prod(n, Q, Jv, x);

%use power method with parameters
eps = 10^(-15);
k_max = 1000;

%for initial condition x = e
e = ones(n,1);
[lambda1, x1] = power_method(matrix_prod, e, eps, k_max);

%for other initial condition
x0 = [1; -1; 1; -1; 1; -1; 1; -1; 1; -1;];
[lambda2, x2] = power_method(matrix_prod, x0, eps, k_max);

%compute relative eigenvalue residual
res1 = max(abs(matrix_prod(x1) - lambda1*x1))/max(abs(x1));
res2 = max(abs(matrix_prod(x2) - lambda2*x2))/max(abs(x2));

%rank websites
[PageRank1, I1] = sort(x1, 'descend');
[PageRank2, I2] = sort(x2, 'descend');
I1(1:10)
I2(1:10)