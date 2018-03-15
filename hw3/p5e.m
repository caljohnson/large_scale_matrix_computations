%Problem 5 - Part E
%MAT226B - HW 3
%Carter Johnson

%use pcg to solve Ax=b for A,b in pcg_small.mat
%using CG without preconditioning 
%and CG with incomplete Chol preconditioning

load('pcg_small.mat');
n = size(A,1);
format long e

[J,I] = find_JI(A)
VL = incomplete_chol(A, J, I)

%initial guess
x0 = ones(n,1);
%pcg tolerance
tol = 10^(-9);

x_uncond = pcg(A,b,tol,[],[],[],x0)
rel_res_norm = norm(b-A*x_uncond)/norm(b)

%make incomplete Chol factors "inverse" functions
M1 = @(x) lower_tri_system_solve(J, I, VL, x);
M2 = @(b) lower_tri_transpose_system_solve( J, I, VL, b );
x_cond = pcg(A,b,tol,[],M1,M2,x0)
rel_res_norm2 = norm(b-A*x_cond)/norm(b)