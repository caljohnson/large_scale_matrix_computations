%Test Rectangular Poisson Solver

%tests to make sure its working
m = 8;
h = 1/(2^m); % grid spacing

ic = 3; % x interval length
id = 1; % y interval length

mx = ic/h - 1; %number of interior x grid points 
my = id/h - 1; %number of interior y grid points

[x,y] = meshgrid(h:h:ic-h, h:h:id-h);

f = exp(-x.^2 - y.^2);

b1 = zeros(1,my);
b0 = ones(1,my);
c1 = zeros(1,mx);
c0 = ones(1,mx);

V = rectangular_twod_poissons(m,f,ic,id,b0,b1,c0,c1);

figure(1); clf;
f = [b0' f b1'];
f = [[0 c0 0]; f; [0 c1 0];];
surf(f); shading interp

figure(2); clf;
V = [b0' V b1'];
V = [[0 c0 0]; V; [0 c1 0];];
surf(V); shading interp
