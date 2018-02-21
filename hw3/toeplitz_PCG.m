function [ x ] = toeplitz_PCG( t, b, tol, x0)
%TOEPLITZ_PCG Solves linear systems C_T^(-1)Tx = C_T^(-1)b
%   where T is a Toeplitz matrix and C_T is the circulant matrix that
%   minimizes ||C - T||_F
% Input: t = [t-(n-1) t-1 ... t0 t1 ... tn-1] - row vector defining Toeplitz matrix T
%        b - column vector
%        tol - tolerance for pcg method
%        x0 - initial guess for pcg method
% Output: x - solution vector to C_T^(-1)Tx = C_T^(-1) b

T_prod = @(y) toeplitz_matrix_product(t, y); %function to multiply T*y, y column vector
% T_prod(b)

%get n
n = size(b,1);

%construct C_T, c_i = (i t[-(n-i)] - (n-i) t[i] ) / n for i=0,...,n-1
%           C = [c0 c1 ... cn-1
%                cn-1 c0 ... cn-2
%                   ...            
%                c1 ... cn-1 c0 ]
%but t[-(n-1)] = t(1) here, so we need to shift by +n
c = zeros(1,n);
%handle i=0 case separately for indexing reasons
c(1,1) = t(n);
for i = 1:n-1
    c(1,i+1) = (i*t(i) + (n-i)*t(n+i) ) / n;
end
%set up as c = [c_0 c_n-1 ... c_2 c_1] - first COLUMN of C_T
c = [c(1,1) c(1,n:-1:2)];

%function to follow uses first COLUMN of C_T
CTinv_prod = @(y) circulant_inverse_matrix_product(c, y); %function to multiply C_T^(-1)y, y col vector


x = pcg(T_prod, b, tol,[], CTinv_prod,[], x0);

end

