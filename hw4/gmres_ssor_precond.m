function [  x, rho, its ] = gmres_ssor_precond( A, b, k0, d_Type)
%GMRES_SSOR_PRECOND Runs GMRES w/ and w/o restarts, tracking residuals
%   with SSOR-type preconditioning to solve Ax=b
%   Uses GMRES function and displays residuals at each restart
%   INPUT: matrix A
%          vector b
%          restart parameter k0
%          d_Type - decides to use D=D0 (==1) or D=10 I (o/w)

n = size(b,1); %get size
x0 = ones(n,1); %set x0
tol = 1e-8; %convergence tolerance
maxits = 100; %maximum outer iterations

%set up matrix-vector product q' = A'v'
D_0 = diag(diag(A));
F = -tril(A,-1);
G = -triu(A,-1);

if d_Type == 1
    D = D_0;
else
    D = 10*speye(size(A));
end

aprime = @(v) matrix_vector_prod_aprime(v,D,D_0, F,G);

if k0 == 0
    k0 = [];
end

%solve M1 b' = b-Ax0 for b'
% (D-F) * D^(-1) b' = b-Ax0
[JL, IL, VL] = coo_to_csc_lower_tri((D-F)*D^(-1)); %convert sparse matrices to CSC
bprime = lower_tri_system_solve(JL, IL, VL, b-A*x0);

%set x0' = 0 and run gmres on A'x'=b' for x'
%set r0' = b'
r0prime = norm(bprime,2);
[xprime,~,~,~,resvec] = gmres(aprime,bprime,k0,tol, maxits);

%solve M2 w = x' for w and set x = x0 + w
%  (D-G)w  = x'
[JU, IU, VU] = coo_to_csc_upper_tri(D-G); %convert sparse matrices to CSC
x = x0 + upper_tri_system_solve(JU, IU, VU, xprime);

%compute relative residual norms
rho = resvec./r0prime;
%compute total # of matrix-vector products q=Av
its = length(rho);

end
