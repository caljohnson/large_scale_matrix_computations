%Problem 3 - HW 4
%MAT 226B
%Carter Johnson

%solve Ax = b to relative residual norm of tol = 1e ? 9
%Run full GMRES, as well as restarted GMRES with restart parameters 
%k0 = 2, 5, 10, 20, 50, 100, all without using preconditioning. 
%Produce a single graph that shows log ?k, k = 0, 1, . . . , for all your runs. 
%get list showing the total number of matrix-vector products for all your runs.

load('HW4_Problem3.mat');

k0 = [0; 2; 5; 10; 20; 50; 100;];

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
legend('full GMRES', 'k_0 = 2','k_0 = 5', 'k_0 = 10', 'k_0 = 20', 'k_0 = 50', 'k_0 = 100');
xlabel('GMRES steps'); ylabel('\rho_k');

table(k0, its_list)