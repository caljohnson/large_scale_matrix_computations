function [ J, I, VU ] = coo_to_csc_upper_tri( U )
%COO_TO_CSC Generates CSC sparse matrix from 
%           an upper-triangular COO sparse matrix U
%   constructs integer vectors J and I and real vector V_L
%   to store U in CSC format

%construct pointer U for upper triangular matrix
[J, K, VU] = find(U);
I = [1; 1+find(J-K==0);];
end

