function [ A1 ] = make_skew_symm( m )
%MAKE_SKEW_SYMM Makes the skew-symmetric matrix A1
%   corresponding to the discretization of
%   u_x  = (u(x+1,y) - u(x-1,y))/2h

h = 1/(m+1);
e = ones(m,1);
A = spdiags([-1.*e 0.*e 1.*e], -1:1, m, m);
A(1,1) = -1;
A(end,end) = 1;
A(1,:) = A(1,:)./h;
A(2:end-1,:) = A(2:end-1,:)./(2*h);
A(end,:) = A(end,:)./h;
                                                 
Ar = repmat(A, 1, m);                                   % Repeat Matrix m times
Ac = mat2cell(Ar, size(A,1), repmat(size(A,2),1,m));    % Create Cell Array Of Orignal Repeated Matrix
A1 = blkdiag(Ac{:});                                    % Desired Result
end

