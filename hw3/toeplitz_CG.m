function [ x ] = toeplitz_CG( t, b, tol, x0 )
%TOEPLITZ_CG Solves linear systems Tx = b
%   where T is a Toeplitz matrix, using CG method
% Input: t = [t-(n-1) t-1 ... t0 t1 ... tn-1] - row vector defining Toeplitz matrix T
%        b - column vector
%        tol - tolerance for pcg method
%        x0 - initial guess for pcg method
% Output: x - solution vector to Tx = b

T_prod = @(y) toeplitz_matrix_product(t, y); %function to multiply T*y, y column vector
x = pcg(T_prod, b, tol, [],[],[], x0);

end

