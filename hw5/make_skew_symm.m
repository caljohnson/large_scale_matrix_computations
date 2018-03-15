function [ A1 ] = make_skew_symm( m )
%MAKE_SKEW_SYMM Makes the skew-symmetric matrix A1
%   corresponding to the discretization of
%   u_x  = (u(x+1,y) - u(x-1,y))/2h

h = 1/(m+1);
e = ones(m^2,1);
A1 = spdiags([-1.*e 0.*e 1.*e], -1:1, m^2, m^2)./(2*h);
end

