function [ V ] = twod_poissons( m, f, b0, b1, c0, c1)
%TWOD_POISSONS 2-D Poisson's Eqn Solver using DFT
%   solves the 2-D Poisson eqn with general BCs using DFT
%   INPUT:  m - discretization size
%           F - matrix of the RHS function at grid points (x_j, y_k)
%           b0, b1, c0, c1 - row-vectors of boundary functions
%   OUTPUT: V - matrix of solution values at grid points

h = 1/(m+1); %mesh size

%add BC's to the RHS function
f(1,:) = f(1,:) + c0/h^2;
f(end,:) = f(end,:) + c1/h^2;
f(:,1) = f(:,1) + b0'/h^2;
f(:,end) = f(:,end) + b1'/h^2;

%compute F' = Z'FZ = Z(ZF)'
f_prime = zeros(m,m);
for j = 1:m
  f_prime(j,:) = fft_matrix_vector_prod(f(:,j),h)';
end
for j = 1:m
  f_prime(:,j) = fft_matrix_vector_prod(f_prime(:,j),h);
end


%compute V_p,jk = h^2 f_jk /(lambda_j + lambda_k)
%eigenvalues
k = 1:m;
lambdas = 2*(1-cos(pi*h*k));
lambda_sums = zeros(m,m);
for j = 1:m
   lambda_sums(j,1:m) = lambdas(j)+lambdas;
end
V_p = h^2*f_prime./lambda_sums;

%compute V = Z V_p Z' = Z (Z V_p)'
V = zeros(m,m);
for j = 1:m
  V(j,:) = fft_matrix_vector_prod(V_p(:,j),h)';
end
for j = 1:m
  V(:,j) = fft_matrix_vector_prod(V(:,j),h);
end
end

