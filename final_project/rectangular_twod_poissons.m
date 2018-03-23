function [ V ] = rectangular_twod_poissons( m, f, ic, id, b0, b1, c0, c1)
%RECTANGULAR_TWOD_POISSONS 2-D Poisson's Eqn Solver using DFT on a
%                                               rectangular domain
%   solves the 2-D Poisson eqn with general BCs using DFT
%   INPUT:  m - discretization factor (size = 2^m)
%           f - matrix of the RHS function at grid points (x_j, y_k) of
%                               size (my,mx)
%           ic - length of x-interval (c, c+ic) defining rectangular domain
%           id - length of y-interval (d, d+id) defining rectangular domain
%           b0, b1, c0, c1 - row-vectors of boundary functions
%   OUTPUT: V - matrix of solution values at grid points

h = 1/(2^m); %get grid spacing
mx = ic/h - 1; %number of x grid points
my = id/h - 1; %number of y grid points
% x = h:h:ic-h;
% y = h:h:id-h;

%add BC's to the RHS function
f(1,:) = f(1,:) + c0/h^2;
f(end,:) = f(end,:) + c1/h^2;
f(:,1) = f(:,1) + b0'/h^2;
f(:,end) = f(:,end) + b1'/h^2;

%compute F' = Z'FZ = Z(ZF)'
f_prime = twod_fft_matrix(f);

%compute V_p,jk = h^2 f_jk /(lambda_j + lambda_k)
% denom = 4/h^2*(sin((x'*pi)/(2*(mx+1))).^2*ones(1,my) + ...
%                         ones(mx,1)*sin((y*pi)/(2*(my+1))).^2);
% V_p = f_prime./denom';

%eigenvalues
kx = 1:mx;
ky = 1:my;
lambdas_x = 2*(1-cos(pi*kx/(mx+1)));
lambdas_y = 2*(1-cos(pi*ky/(my+1)));
lambda_sums = zeros(my,mx);
if mx <= my
    for j = 1:mx
        lambda_sums(1:my, j) = lambdas_x(j)+lambdas_y;
    end
else
    for j = 1:my
        lambda_sums(j, 1:mx) = lambdas_x+lambdas_y(j);
    end
end
V_p = h^2*f_prime./lambda_sums;

%compute V = Z V_p Z' = Z (Z V_p)'
V = twod_fft_matrix(V_p);

end

