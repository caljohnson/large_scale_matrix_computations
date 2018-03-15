function [ x ] = lower_tri_transpose_system_solve( J, I, VL, b )
%LOWER_TRI_TRANSPOSE_SYSTEM_SOLVE Solves L'x = b
%   where L is given in CSC format by J, I, and VL
%       and is lower-triangular
% Input: vector J of row indices
%        vector I of column pointers
%        vector VL of sparse values
%        vector b
% Output: x - soln to L'x = b

%preallocate space for x
n = size(b,1);
% x = zeros(n,1);

%backsolve row-wise, since L' is given row-wise in J,I, VL
for j=n:-1:2
    b(j) = b(j)/VL(I(j));
    %find VL indices of stuff in ROW j of L
    inds = find(J==j);
    %for each ind do
    for i = 1:length(inds)-1
       %find corresponding column # of L
       col_no = find(I - inds(i)<=0,1,'last');
       b(col_no) = b(col_no) - b(j)*VL(inds(i));
    end
end
b(1) = b(1)/VL(1);
x=b;
end

