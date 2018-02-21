function [ y ] = circulant_inverse_matrix_product( c , x )
%CIRCULANT_INVERSE_MATRIX_PRODUCT Matrix-vector product C^(-1) x
%   where C is a circulant matrix in R^(n x n):
%           C = [c0 cn-1 ... c1
%                c1 c0 ... cn-1
%                   ...            
%                cn-1 .. c1 c0 ]
%       and x is a vector in R^n
%Inputs:  c = [c0 c1 ... cn-1] - row vector defining circulant matrix C
%                                  (c is the first column of C)
%         x  - vector to mutliply by
%Outputs: y = Cx

d = fft(c);  %obtain eigenvalues of C
d = 1./d; %set to eigenvalues of C^(-1), which are 1/eigvals of C
D = spdiags(d(:),0,size(d,2),size(d,2));

y = ifft(full(D*fft(x))); %y = C^(-1)x = (F'DF)x since C^(-1) = F'DF
end

