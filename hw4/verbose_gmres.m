function [ x, rho, its] = verbose_gmres( A, b, k0 )
%VERBOSE_GMRES Runs GMRES w/ and w/o restarts, tracking residuals
%   Uses GMRES function and displays residuals at each restart

n = size(b,1); %get size
x0 = ones(n,1); %set x0
tol = 1e-8; %convergence tolerance
maxits = 500; %maximum outer iterations

r0 = norm(A*x0 - b,2);

if k0 == 0
    k0 = [];
end
[x,~,~,~,resvec] = gmres(A,b,k0,tol, maxits);

%compute relative residual norms
rho = resvec./r0;
%compute total # of matrix-vector products q=Av
its = length(rho);

end

