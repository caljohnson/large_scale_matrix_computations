function [ J, I ] = find_JI( A )
%FIND_JI Computes integer vectors J, I from sparse symmetric A
%  INPUT: Sparse symmetric matrix A
%  OUTPUT: Row-indices vector J, Pointer vector I
[J,K] = find(tril(A));
I = find(J-K==0);
I = [I; nnz(tril(A)) + 1;];
end

