function [ c ] = lower_tri_system_solve( J, I, VL, b )
%LOWER_TRI_SYSTEM_SOLVE Solves Lc = b
%   where L is a lower triangular matrix given in CSC format by J, I, and VL
% Input: vector J of row indices
%        vector I of column pointers
%        vector VL of sparse values
%        vector b
% Output: c - soln to Lc = b

%allocate store for c
n = size(b,1);
c = zeros(n,1);
%forward solve by columns
for k=1:n-1;
    c(k) = b(k)/VL(I(k));
    b(J(I(k):I(k+1)-1)) = b(J(I(k):I(k+1)-1)) - VL(I(k):I(k+1)-1).*c(k); 
end
c(end) = b(end)/VL(I(n));

end

