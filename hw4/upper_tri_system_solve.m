function [ x ] = upper_tri_system_solve( J, I, VU, c )
%UPPER_TRI_SYSTEM_SOLVE Solves Ux = c
%   where U is given in CSC format by J, I, and VU
% Input: vector J of row indices
%        vector I of column pointers
%        vector VU of sparse values
%        vector c
% Output: x - soln to Ux = c

%preallocate space for x
n = size(c,1);
x = zeros(n,1);

%backsolve columnwise
for k=n:-1:2
    x(k) = c(k)/VU(I(k+1)-1);
    J_indeces = J(I(k):I(k+1)-2); 
    c(J_indeces) = c(J_indeces) - VU(I(k):I(k+1)-2).*x(k); 
end
x(1) = c(1)/VU(I(2)-1);

end
