%Part 3 - C
%MAT226B  - Final Project
%Carter Johnson


m = 3;
h = 1/(2^m);

%case 4
alpha = 2;
beta = 7/2;
gamma = 2;
a = 0.25;
b = 1;

f= @(x,y) sin(beta*pi*x).*(beta^2*pi^2*y.^alpha.*cos(gamma*pi*y) + ...
    gamma^2*pi^2*y.^alpha.*cos(gamma*pi*y) - ...
    alpha*(alpha-1)*y.^(alpha-2).*cos(gamma*pi*y) + ...
    2*alpha*gamma*pi*y.^(alpha-1).*sin(gamma*pi*y)); 

g = @(x,y) y.^alpha.*sin(beta*pi.*x).*cos(gamma*pi.*y);

u = @(x,y) y.^alpha.*sin(beta*pi.*x).*cos(gamma*pi.*y);
[x1,y1] = meshgrid(h:h:1-h,h:h:3-h); 
[x2,y2] = meshgrid(a+h:h:a+4-h,b+h:h:b+1-h); 
u1 = u(x1,y1);
u2 = u(x2,y2);
figure(1); clf;
surf(x1,y1,u1); shading interp; hold on;
surf(x2,y2,u2); shading interp; hold off;


[A,c, I1, I2 ] = generate_data( f,g, a,b, m);
v0_1 = repermute_f_R1(u1',h,a,b);
v0_2 = repermute_f_R2(u2',h,a);
v0 = [v0_1(1:end-((1-a)/h-1)*(1/h-1)); v0_2];
tol = 1e-10;
nmax = 1000;

% v = multiplicative_schwarz_method(A,c,I1,I2,m, a,b,v0,tol,nmax);
v = additive_schwarz_pcg(A,c,I1,I2,m, a,b,v0,tol,nmax);

v1 = permute_f_R1(I1*v,h,a,b,1,3);
v2 = permute_f_R2(I2*v,h,a,4,1);

max1 = max(max(abs(v1-u1')))/max(max(abs(u1)));
max2 = max(max(abs(v2-u2')))/max(max(abs(u2)));
max(max1,max2)

display('done');

figure(2); clf;
surf(x1,y1,v1'); shading interp; hold on;
surf(x2,y2,v2'); shading interp; hold off;
% title('Computed'); zlabel('V'); xlabel('x'); ylabel('y');
