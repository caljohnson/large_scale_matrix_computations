function [v, RELRES, its ] = generalized_2d_poissons( gamma, m, f, g)
%GENERALIZED_2D_POISSONS Solves generalized 2D Poissons eqn
%       -uxx - uyy + gamma*ux = f(x, y), (x, y) ? R := (0, 1) × (0, 1),
%           u = g(x, y), (x, y) ? ?R.
%INPUT: gamma, m - discretization size, f - RHS function, g - boundary function
%OUTPUT: approximate solution v, relative residual norm RELRES

%create meshgrid on (0,1)x(0,1) with discretization size m
h = 1/(m+1);
[x,y] = meshgrid(h:h:1-h, h:h:1-h);

%sample function f on meshgrid points to get discrete RHS matrix F
F = f(x,y);

%sample function g on boundaries to set boundary vectors b0, b1, c0, c1
b0 = g(x(1,:), 0);
b1 = g(x(1,:), 1);
c0 = g(0, y(1,:));
c1 = g(1,y(1,:));

%create skew-symm matrix
A1 = make_skew_symm(m);
%create matrix-vector product A'x = (I + gamma*M_1^(-1) A_1)x 
A_prime = @(x) x + gamma*reshape(twod_poissons(m, reshape(A1*x, [m,m]), ...
                zeros(1,m),zeros(1,m), zeros(1,m), zeros(1,m)), [m^2,1]);

%hit RHS with left-preconditioner
b_prime = reshape(twod_poissons(m, F, b0, b1, c0, c1), [m^2,1]);

%run gmres
[v,~,RELRES, its] = gmres(A_prime,b_prime, [], 1e-10, 1000);

%return total iteration count
its = its(2);
end

