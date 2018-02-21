function [ y ] = circulant_matrix_product( c , x )
%CIRCULANT_MATRIX_PRODUCT Matrix-vector product C x
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
D = spdiags(d(:),0,size(d,2),size(d,2)); %make sparse diagonal matrix w/ eigenvals of C on diag

y = ifft(D*fft(x)); %y = Cx = (F'DF)x since C = F'DF
end

