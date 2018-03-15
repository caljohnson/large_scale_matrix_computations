%Problem 5 - HW 4
%MAT 226B
%Carter Johnson

%solve Ax = b to relative residual norm of tol = 1e ? 9
%Run full GMRES, as well as restarted GMRES with restart parameters 
%k0 = 5, 10, 20
%for no-preconditioning, diagonal right-preconditioning, and SSOR-type
%preconditioning w/ D=D_0 and D=10 I
%For each run, produe a single graph that shows log rho k, k = 0, 1, . . . . 
%and get list showing the total number of matrix-vector products for all your runs.

% load('HW4_Problem5b_1.mat');
load('HW4_Problem5b_2.mat');

k0 = [0; 5; 10; 20;];

%store rel. residual vectors
rhos = [];
%store list of total m-v-prods number
its_list = zeros(size(k0));
%store list of solns xk
xk_list = zeros(size(k0,1), size(b,1));

figure(1); clf;
for ii =1:size(k0,1)
   [xk, rho, its] = verbose_gmres(A,b,k0(ii));
   
   its_list(ii) = its;
   xk_list(ii, :) = xk';
   plot(log(rho)); hold on;
end
hold off;
title('No Preconditioning');
legend('full GMRES', 'k_0 = 5', 'k_0 = 10', 'k_0 = 20');
xlabel('GMRES steps'); ylabel('\rho_k');
table(k0, its_list)

figure(2); clf;
for ii =1:size(k0,1)
   [xk, rho, its] = gmres_diag_precond(A,b,k0(ii));
   
   its_list(ii) = its;
   xk_list(ii, :) = xk';
   plot(log(rho)); hold on;
end
hold off;
title('Diagonal Right-Preconditioning');
legend('full GMRES', 'k_0 = 5', 'k_0 = 10', 'k_0 = 20');
xlabel('GMRES steps'); ylabel('\rho_k');
table(k0, its_list)

figure(3); clf;
for ii =1:size(k0,1)
   [xk, rho, its] = gmres_ssor_precond(A,b,k0(ii),1);
   
   its_list(ii) = its;
   xk_list(ii, :) = xk';
   plot(log(rho)); hold on;
end
hold off;
title('SSOR-Type Preconditioning, D=D_0');
legend('full GMRES', 'k_0 = 5', 'k_0 = 10', 'k_0 = 20');
xlabel('GMRES steps'); ylabel('\rho_k');
table(k0, its_list)

figure(4); clf;
for ii =1:size(k0,1)
   [xk, rho, its] = gmres_ssor_precond(A,b,k0(ii),0);
   
   its_list(ii) = its;
   xk_list(ii, :) = xk';
   plot(log(rho)); hold on;
end
hold off;
title('SSOR-Type Preconditioning, D=10 I');
legend('full GMRES', 'k_0 = 5', 'k_0 = 10', 'k_0 = 20');
xlabel('GMRES steps'); ylabel('\rho_k');
table(k0, its_list)