%Problem 1 - HW 5
%MAT 226B
%Carter Johnson

%Test generalized 2-D Poisson's eqn solver
%       ?uxx ? uyy + ?ux = f(x, y), (x, y) ? R := (0, 1) × (0, 1),
%           u = g(x, y), (x, y) ? ?R.
% using left-preconditioner M_1 = A_0 = discrete 2D Laplacian
% and GMRES w/ tol 1e-10

%test RHS function
f = @(x,y) x.^3 .* y.^2 .* exp(2-x-y);
%test boundary function
g = @(x,y) 0.*x + 0.*y + 1;

%gamma values
gamma = [1; 10; 50; 100; 1000;];

%single m value
m=100;

RELRES = zeros(5,1);
its = zeros(5,1);
for jj = 1:size(gamma,1);
    %run generalized 2d poisson solver
%     gamma(jj)
    [v, RELRES(jj), its(jj)] = generalized_2d_poissons(gamma(jj),m,f,g);
end
plot(gamma, its);
table(gamma, RELRES, its);
