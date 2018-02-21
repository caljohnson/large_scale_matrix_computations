%Problem 2 - test cases
%MAT226B  - HW 3
%Carter Johnson

a = 1;
b = 2;
mx = 27;
my = 2*27;
hx = a/(mx+1);
hy = b/(my+1);
[x,y] = meshgrid(hx:hx:a-hx, hy:hy:b-hy);

%case 1
alpha = 0;
beta = 1;
gamma = 0.5;
f= beta^2*pi^2*sin(beta*pi*x).*cos(gamma*pi*y) + gamma^2*pi^2*sin(beta*pi*x).*cos(gamma*pi*y);
f = f';
b0 = sin(beta*pi*(hx:hx:a-hx));
b1 = zeros(1,mx);
c0 = zeros(1,my);
c1 = zeros(1,my);

% %case 2
% alpha = 1;
% beta = 1.5;
% gamma = 2;
% f= (beta^2 + gamma^2)*pi^2*y.*sin(beta*pi*x).*cos(gamma*pi*y) + (2*gamma*pi)*sin(beta*pi*x).*sin(gamma*pi*y);
% f = f';

% %case 3
% alpha = 2;
% beta = 3;
% gamma = 0.5;

% %case 4
% alpha = 5;
% beta = 3;
% gamma = 1;

% %case 5
% alpha = 5;
% beta = 5;
% gamma = 3;

% %cases 3-5
% f= sin(beta*pi*x).*(beta^2*pi^2*y.^alpha.*cos(gamma*pi*y) + ...
%     gamma^2*pi^2*y.^alpha.*cos(gamma*pi*y) - ...
%     alpha*(alpha-1)*y.^(alpha-2).*cos(gamma*pi*y) + ...
%     2*alpha*gamma*pi*y.^(alpha-1).*sin(gamma*pi*y)); 
% f = f';
figure(1); clf;
subplot(1,3,1);
surf(f); shading interp;
title('RHS'); zlabel('f'); xlabel('x'); ylabel('y');
% 
% %cases 2-5
% b0 = 0^alpha*sin(beta*pi*(h:h:1-h))*cos(gamma*pi*0);
% b1 = 1^alpha*sin(beta*pi*(h:h:1-h))*cos(gamma*pi*1);
% c0 = (h:h:1-h).^alpha*sin(beta*pi*0).*cos(gamma*pi*(h:h:1-h));
% c1 = (h:h:1-h).^alpha*sin(beta*pi*1).*cos(gamma*pi*(h:h:1-h));

V_actual = y.^alpha.*sin(beta*pi*x).*cos(gamma*pi*y);
V_actual = V_actual';
subplot(1,3,2);
surf(V_actual); shading interp;
title('Analytic'); zlabel('V'); xlabel('x'); ylabel('y');

V_test = generalized_twod_poissons(mx, my, a,b,f,b0,b1,c0,c1);
subplot(1,3,3);
surf(V_test); shading interp;
title('Computed'); zlabel('V'); xlabel('x'); ylabel('y');

max(max(abs(V_test - V_actual)))
m
