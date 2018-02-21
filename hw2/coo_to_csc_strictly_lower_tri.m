function [ J, I, VL ] = coo_to_csc_strictly_lower_tri( L )
%COO_TO_CSC Generates CSC sparse matrix from 
%           a unit lower-triangular COO sparse matrix L
%   constructs integer vectors J and I and real vector V_L
%   to store the strictly lower-triangular part of L in CSC format

sl_L = tril(L,-1);

%construct pointer I for lower triangular matrix
[J, K, V] = find(L);
I = find(J-K==0);

%convert to strictly lower triangular
%by removing unit diagonal elements from J, I, VL
VL = V(J-K~=0); %remove corresponding diagonal values from V
J = J(J-K~=0); %remove all diagonal entries from J
I(2:end) = I(2:end)-(1:size(I,1)-1)'; %subtract 1 from each nnz count in I for each additional row
I = [I; nnz(sl_L)+1];
end

