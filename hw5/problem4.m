%Problem 4 - HW 5
%MAT 226B
%Carter Johnson

%Employ your implementation of the Arnoldi process
%to compute approximate eigenpairs of the matrix
%A = A(gamma) = A' = I + gamma A0^{-1}A 

%gamma values
gamma = [1; 10; 50; 100; 1000;];

m=100;
%given kmax
kmax = 300;

%given r
load('HW5_P4.mat');

%set up storage for min/max rhos
min_rhos = zeros(size(gamma));
min_rhos_eig = zeros(size(gamma));
max_rhos = zeros(size(gamma));
max_rhos_eig = zeros(size(gamma));

for jj = 1:size(gamma,1);
    %create skew-symm matrix
    A1 = make_skew_symm(m);
   
    %create matrix-vector product A'x = (I + gamma*M_1^(-1) A_1)x 
    A_prime = @(x) x + gamma(jj)*reshape(twod_poissons(m, reshape(A1*x, [m,m]), ...
                zeros(1,m),zeros(1,m), zeros(1,m), zeros(1,m)), [m^2,1]);
    %precondition residual
%     r = reshape(twod_poissons(m, reshape(r, [m,m]),zeros(1,m), zeros(1,m),zeros(1,m), zeros(1,m)), [m^2,1]);
    
    %compute Arnoldi vectors and Hk
    [V,H] = arnoldi_process(A_prime, r, kmax);
    
    %compute eigenpairs of Hk
    [Z, D] = eig(H(1:kmax, 1:kmax));
    
    %compute residual norms of approx eigenpairs
    rhos = H(kmax+1, kmax)*abs(Z(end,:));
    
    %print min rho and max rho
    [~,ii] = min(rhos);
    [~, II] = max(rhos);
    min_rhos(jj) = rhos(ii);
    min_rhos_eig(jj) = D(ii,ii);
    max_rhos(jj) = rhos(II);
    max_rhos_eig(jj) = D(II,II);
    %plot all approx eigenvals
    figure(jj); clf;
    plot(diag(D), 'o'); 
    pltitle = strcat('\gamma = ', num2str(gamma(jj)));
    xlabel('Re \lambda'); ylabel('Im \lambda'); title(pltitle);
end
format long e
table(gamma, min_rhos, min_rhos_eig)
table(gamma, max_rhos, max_rhos_eig)