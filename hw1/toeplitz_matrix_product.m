function [ y ] = toeplitz_matrix_product( t , x )
%TOEPLITZ_MATRIX_PRODUCT Matrix-vector product T x
%   where T is a Toeplitz matrix in R^(n x n):
%           T = [t0 t1... tn-1
%                t-1 t0 ... tn-1
%                   ...          t1  
%                t-(n-1) .. t-1 t0 ]
%       and x is a vector in R^n
%   by converting creating a circulant matrix C in R^(2n-1 x 2n-1) s.t. T is the leading
%   principal nxn submatrix of C
%Inputs:  t = [t-(n-1) t-1 ... t0 t1 ... tn-1] - row vector defining Toeplitz matrix T
%         x  - (n times 1 ) vector to mutliply by
%Outputs: y = Tx

n = (size(t,2)+1)/2; %get dimension of space R^n: size(t) = 2n-1
x = [x; zeros(n-1,1);]; %add extra zeros to x so we only multiply by T

%c = [t0 t-1 ... t-(n-1) tn-1 ... t1]
c = zeros(1, 2*n-1); %set up 1x2n-1 row vector to define C
c(1,1) = t(n); %t(n) = t_0 -> first element in c
c(1,2:n) = t(n-1:-1:1); %next elements in c = [t-1, ..., t-(n-1)]
c(1,n+1:end) = t(end:-1:n+1); %next elements in c = [tn-1, ..., t1]

y = circulant_matrix_product(c, x);
y = y(1:n);
end

