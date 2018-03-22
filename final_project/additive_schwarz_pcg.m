function [ v ] = additive_schwarz_pcg( A, c, I1,I2,m, a,b, v0,tol,nmax)
%ADDITIVE_SCHWARZ_PCG Uses PCG with Additive Schwarz as preconditioner
%   Input: Discretization A,c
%          Restriction operators I1:R->R1, I1:R->R2
%          Initial guess v0, tolerance tol, its max nmax
%          system info m, a,b

h=1/(2^m); %gridsize

%B1 = I1'*A1^(-1)*I1;
ic = 1;
id = 3;
B1 = @(w) I1'*repermute_f_R1(rectangular_twod_poissons(m, ...
       permute_f_R1(I1*w, h,a,b,ic,id)',ic,id, zeros(1,id/h-1), ...
       zeros(1,id/h-1),zeros(1,ic/h-1),zeros(1,ic/h-1))', h, a,b);

%B2 = I2'*A2^(-1)*I2;
ic = 4;
id = 1; 
B2 = @(w) I2'*repermute_f_R2(rectangular_twod_poissons(m, ...
       permute_f_R2(I2*w, h,a,ic,id)',ic,id, zeros(1,id/h-1), ...
       zeros(1,id/h-1),zeros(1,ic/h-1),zeros(1,ic/h-1))', h, a);
   
%use M^(-1) = B1+B2 as left-preconditioner => A'= M^-1Av = (B1+B2)Av
Aprime = @(w) B1(A*w) + B2(A*w);
%precondition RHS - M^-1c = (B1+B2)c
bprime = B1(c) + B2(c);

%return final PCG iterate
v = pcg(Aprime, bprime, tol, nmax, [], [], v0);

end

