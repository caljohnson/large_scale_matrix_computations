%Test Rectangular Poisson Solver

%tests to make sure its working
m = 6;
h = 1/(2^m); % grid spacing

ic = 3; % x interval length
id = 1; % y interval length

mx = ic/h - 1; %number of interior x grid points 
my = id/h - 1; %number of interior y grid points

alpha = 5;
beta = 3;
gamma = 1;

f= @(x,y) sin(beta*pi*x).*(beta^2*pi^2*y.^alpha.*cos(gamma*pi*y) + ...
    gamma^2*pi^2*y.^alpha.*cos(gamma*pi*y) - ...
    alpha*(alpha-1)*y.^(alpha-2).*cos(gamma*pi*y) + ...
    2*alpha*gamma*pi*y.^(alpha-1).*sin(gamma*pi*y)); 

g = @(x,y) y.^alpha.*sin(beta*pi.*x).*cos(gamma*pi.*y);

u = @(x,y) y.^alpha.*sin(beta*pi.*x).*cos(gamma*pi.*y);

[x,y] = meshgrid(h:h:ic-h, h:h:id-h);
b1 = g(1,h:h:id-h);
b0 = g(0,h:h:id-h);
c1 = g(h:h:ic-h,1);
c0 = g(h:h:ic-h,0);

F = f(x,y);

V = rectangular_twod_poissons(m,F,ic,id,b0,b1,c0,c1);
real_v = u(x,y);

max(max(abs(V - real_v)))

figure(1); clf;
surf(x,y,real_v); shading interp;
figure(2); clf;
surf(x,y,V); shading interp; 

% figure(1); clf;
% f = [b0' f b1'];
% f = [[0 c0 0]; f; [0 c1 0];];
% % surf(f); shading interp
% 
% figure(2); clf;
% V = [b0' V b1'];
% V = [[0 c0 0]; V; [0 c1 0];];
% surf(V); shading interp
