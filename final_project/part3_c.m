%Part 3 - C
%MAT226B  - Final Project
%Carter Johnson

m = 5;
h = 1/(2^m);

schwarz_method = 2; %1 for Multiplicative, 2 for Additive PCG

% %case 1
% case_no = 1;
% a=0.5;
% b=0;
% alpha=0;
% beta=1;
% gamma=0.5;

% %case 2
% case_no = 2;
% a = 0.5;
% b = 0;
% alpha = 2;
% beta = 7/2;
% gamma = 2;

% %case 3
% case_no = 3;
% a=0.25;
% b=1;
% alpha=0;
% beta=1;
% gamma=0.5;

% %case 4
% case_no = 4;
% a = 0.25;
% b = 1;
% alpha = 2;
% beta = 7/2;
% gamma = 2;

% %case 5
% case_no = 5;
% a=0.5;
% b=2;
% alpha=0;
% beta=1;
% gamma=0.5;

%case 6
case_no = 6;
a=0.5;
b=2;
alpha=2;
beta=7/2;
gamma=2;

%even cases
if mod(case_no,2)==0
    f= @(x,y) sin(beta*pi*x).*(beta^2*pi^2*y.^alpha.*cos(gamma*pi*y) + ...
        gamma^2*pi^2*y.^alpha.*cos(gamma*pi*y) - ...
        alpha*(alpha-1)*y.^(alpha-2).*cos(gamma*pi*y) + ...
        2*alpha*gamma*pi*y.^(alpha-1).*sin(gamma*pi*y));
else
    f=@(x,y) (beta^2 + gamma^2)*pi^2*sin(beta*pi*x).*cos(gamma*pi*y);
end

g = @(x,y) y.^alpha.*sin(beta*pi.*x).*cos(gamma*pi.*y);

%set up grid
[x1,y1] = meshgrid(h:h:1-h,h:h:3-h); 
[x2,y2] = meshgrid(a+h:h:a+4-h,b+h:h:b+1-h);
%sample analytic solution
u1 = g(x1,y1);
u2 = g(x2,y2);
figure(1); clf;
surf(x1,y1,u1); shading interp; hold on;
surf(x2,y2,u2); shading interp; hold off;
title('Sampled, true'); zlabel('V'); xlabel('x'); ylabel('y');


[A,c, I1, I2 ] = generate_data( f,g, a,b, m);
v0_1 = repermute_f_R1(u1',h,a,b);
v0_2 = repermute_f_R2(u2',h,a);

v0 = [v0_1(1:end-((1-a)/h-1)*(1/h-1)); v0_2];

tol = 5e-4;
nmax = 1000;

if schwarz_method == 1
    v = multiplicative_schwarz_method(A,c,I1,I2,m, a,b,v0,tol,nmax);
else
    v = additive_schwarz_pcg(A,c,I1,I2,m, a,b,v0,tol,nmax);
end

v1 = permute_f_R1(I1*v,h,a,b,1,3);
v2 = permute_f_R2(I2*v,h,a,4,1);

max1 = max(max(abs(v1-u1')))/max(max(abs(v1)));
max2 = max(max(abs(v2-u2')))/max(max(abs(v2)));
format long e
max(max1,max2)

display('done');

fig2 = figure(2); clf;
surf(x1,y1,v1'); shading interp; hold on;
surf(x2,y2,v2'); shading interp; hold off;
zlabel('V'); xlabel('x'); ylabel('y');
if schwarz_method ==1
    title('Multiplicative Schwarz');
    figname = strcat(strcat('~/Documents/my_classes/winter_quarter_2018/numerics/writeups/final_project/mult_schwarz_case',num2str(case_no)),'.jpg');
    saveas(fig2,figname)
else
    title('Additive Schwarz PCG'); 
    figname = strcat(strcat('~/Documents/my_classes/winter_quarter_2018/numerics/writeups/final_project/add_schwarz_pcg_case',num2str(case_no)),'.jpg');
    saveas(fig2,figname)
end
