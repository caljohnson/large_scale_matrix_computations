%Problem 1 - HW 2 - MAT226B
%Carter Johnson

%Part a

%Compute the Cholesky factor L of A = 2D laplacian
%for the following 5 cases:
% (i) No reordering of the rows and columns of A;
% (ii) Reordering with symamd;
% (iii) Reordering with colamd;
% (iv) Reordering with symrcm;
% (v) Reordering with colperm.
%For each run, report the number of nonzero entries of L
%and submit a plot of the sparsity structure of L

m0 = 67;
m= m0;
A = make_2d_laplacian(m);
% 
% %(i) No reordering
% L1 = chol(A, 'lower');
% display('Number of nonzeros:');
% nnz(L1)
% figure(1);
% spy(L1);
% 
% % (ii) Reordering with symamd;
% p = symamd(A);
% L2 = chol(A(p,p), 'lower');
% display('Number of nonzeros:');
% nnz(L2)
% figure(2);
% spy(L2);
% 
% % (iii) Reordering with colamd;
% p = colamd(A);
% L3 = chol(A(p,p), 'lower');
% display('Number of nonzeros:');
% nnz(L3)
% figure(3);
% spy(L3);
% 
% 
% % (iv) Reordering with symrcm;
% p = symrcm(A);
% L4 = chol(A(p,p), 'lower');
% display('Number of nonzeros:');
% nnz(L4)
% figure(4);
% spy(L4);
% 
% 
% % (v) Reordering with colperm.
% p = colperm(A);
% L5 = chol(A(p,p), 'lower');
% display('Number of nonzeros:');
% nnz(L5)
% figure(5);
% spy(L5);

for i=0:10
   tic
   m = 2^i * m0;
   A = make_2d_laplacian(m);
   p = symamd(A);
   L2 = chol(A(p,p), 'lower');
   i
   m
   toc
end
