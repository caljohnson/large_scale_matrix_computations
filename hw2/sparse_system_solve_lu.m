function [ x ] = sparse_system_solve_lu( A, b, p0, scaling )
%SPARSE_SYSTEM_SOLVE_LU Solves linear system Ax=b
%   where A is a sparse matrix
%   via sparse LU factorization and sparse back/forward solves
% Input: sparse matrix A
%        vector b
%        column permutation vector p0
%        scaling flag - 1=use scaling, 0=no scaling
% Output: x - solution to Ax=b


tic
b = b(p0);
if scaling==1
    [L, U,p,q,D] = lu(A(:,p0), 'vector');
else
    [L, U,p,q] = lu(A(:,p0), 'vector');
end
toc
n = size(p,2);
tic
[JL, IL, VL] = coo_to_csc_strictly_lower_tri(L);
I = speye(n,n);
P = I(p,:);

if scaling==1
    c = unit_lower_tri_system_solve(JL, IL, VL, P*D^(-1)*b);
else
    c = unit_lower_tri_system_solve(JL, IL, VL, P*b);
end

[JU, IU, VU] = coo_to_csc_upper_tri(U);
d = upper_tri_system_solve(JU, IU, VU, c);
qi(q)= 1:n;
x = d(qi);
toc

%print it all out
format long e
nnz(L)
nnz(U)
norm(b-A(:,p0)*x)/norm(b)
x(10), x(100), x(1000), x(100000), x(200000)

end

