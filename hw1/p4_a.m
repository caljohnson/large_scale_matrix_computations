%HW 1 Problem 4a  - MAT226B WQ1028
%Carter Johnson

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
A(A<=1/10)=0;
Q = sparse(A);
e = ones(10,1);
Jv = [3 6 7];

format long e
y = row_stochastic_sparse_matrix_prod(10, Q, Jv, e)
y_half = row_stochastic_sparse_altered_matrix_prod(0.5, 10, Q, Jv, e)
y_eightfive = row_stochastic_sparse_altered_matrix_prod(0.85, 10, Q, Jv, e)


%second part
load('www0.mat')
load('x0.mat')

y = row_stochastic_sparse_matrix_prod(size(x0,1), Q, Jv, x0);
format long e
display('y')
y(2)
y(222222)
y(300000)
y(400000)

display('y_0.85')
y = row_stochastic_sparse_altered_matrix_prod(0.85, size(x0,1), Q, Jv, x0);
y(2)
y(222222)
y(300000)
y(400000)