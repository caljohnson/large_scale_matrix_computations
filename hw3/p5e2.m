%Problem 5 - Part E part 2
%MAT226B - HW 3
%Carter Johnson

%use pcg to solve Ax=b for A,b in pcg_large.mat
%using CG without preconditioning 
%and CG with incomplete Chol preconditioning

load('pcg_large.mat');
n = size(A,1);
format long e

[J,I] = find_JI(A);
VL = incomplete_chol(A, J, I);

%initial guess
x0 = ones(n,1);
%pcg tolerance
tol = 10^(-9);

% x_uncond = pcg(A,b,tol,[],[],[],x0);
% rel_res_norm = norm(b-A*x_uncond)/norm(b)
% x_uncond(1), x_uncond(10000), x_uncond(100000), x_uncond(200000), x_uncond(262144)

%make incomplete Chol factors "inverse" functions
M1 = @(x) lower_tri_system_solve(J, I, VL, x);
M2 = @(b) lower_tri_transpose_system_solve( J, I, VL, b );
x_cond = pcg(A,b,tol,[],M1,M2,x0);
rel_res_norm2 = norm(b-A*x_cond)/norm(b)

table(x_cond(1), x_cond(10000), x_cond(100000), x_cond(200000), x_cond(262144))
% table(J(1), J(10), J(100), J(1000), J(end))
% table(I(1), I(10), I(100), I(1000), I(end))
% table(VL(1), VL(10), VL(100), VL(1000), VL(end))