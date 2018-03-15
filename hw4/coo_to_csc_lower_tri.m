function [ J, I, VL ] = coo_to_csc_lower_tri( L )
%COO_TO_CSC_LOWER_TRI Computes integer vectors J, I from lower triangular matrix A
%  INPUT: Lower triangular matrix L
%  OUTPUT: Row-indices vector J, Pointer vector I, values VL
[J,K, VL] = find(tril(L));
I = find(J-K==0);
I = [I; nnz(tril(L)) + 1;];
end

