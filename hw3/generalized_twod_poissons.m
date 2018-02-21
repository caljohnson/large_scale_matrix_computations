function [ V ] = twod_poissons( mx,my, a,b, f, b0, b1, c0, c1)
%TWOD_POISSONS 2-D Poisson's Eqn Solver using DFT
%   solves the 2-D Poisson eqn with general BCs using DFT
%   INPUT:  mx - discretization size in x
%           my - discretization size in y
%           a - grid endpoint in x - assumed x in (0,a)
%           b - grid endpoint in y - assumed y in (0,b)
%           F - matrix of the RHS function at grid points (x_j, y_k)
%           b0, b1, c0, c1 - row-vectors of boundary functions
%   OUTPUT: V - matrix of solution values at grid points

hx = a/(mx+1); %mesh size in x
hy = b/(my+1); %mesh size in x

%add BC's to the RHS function
f(1,:) = f(1,:) + c0;
f(mx,:) = f(mx,:) + c1;
f(:,1) = f(:,1) + b0'*(hx^2/hy^2);
f(:,my) = f(:,my) + b1'*(hx^2/hy^2);

%compute F' = Z'FZ = Z(ZF)'
f_prime = zeros(mx,my);
for j = 1:my
  f_prime(:,j) = fft_matrix_vector_prod(f(:,j),hy);
end
for j = 1:mx
  f_prime(j,:) = fft_matrix_vector_prod(f_prime(j,:)',hx)';
end


%compute V_p,jk = h^2 f_jk /(lambda_j + lambda_k)
%eigenvalues
k = 1:mx;
lambdas = 2*(1-cos(pi*hx*k));
lambda_sums = zeros(mx,my);
for j = 1:mx
   lambda_sums(j,1:my) = lambdas(j)+lambdas;
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

