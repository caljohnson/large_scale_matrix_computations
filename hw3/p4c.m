%Problem 4 - part C
%MAT 226B - HW 3
%Carter Johnson

%test on Toeplitz matrices with t_j = 1/(1+sqrt(j))^p 
t_entry = @(j,p) 1./(1+sqrt(j)).^p;

n = 10^6;
% n = 10;

%Cases p
p = [1; 0.1; 0.01;];

%RHS vector
b = ones(n,1);
%initial guess
x0 = zeros(n,1);
%tolerance for pcg
tol = 10^(-9);

x_cg = zeros(n,3);
x_pcg = zeros(n,3);

for i = 1:2
% for i = 1:3
    j = 0:1:n-1;
    t = t_entry(j,p(i));
%     T = toeplitz(t);
    t = [t(n:-1:2) t(1:n)];
%     x = T\b;
    
    x_cg(:,i) = toeplitz_CG(t,b,tol,x0);
%     norm(abs(x - x_cg(:,i)))
    x_pcg(:, i) = toeplitz_PCG(t,b,tol,x0);
%     norm(abs(x - x_pcg(:,i)))
end
format long e

% %case n=10
% table(x_cg)
% table(x_pcg)

%case n = 10^6
disp('CG');
x = x_cg(:,1);
table(x(1), x(100000), x(500000), x(700000), x(1000000))
x = x_cg(:,2);
table(x(1), x(100000), x(500000), x(700000), x(1000000))

disp('PCG');
x = x_pcg(:,1);
table(x(1), x(100000), x(500000), x(700000), x(1000000))
x = x_pcg(:,2);
table(x(1), x(100000), x(500000), x(700000), x(1000000))
