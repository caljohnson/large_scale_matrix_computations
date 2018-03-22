function [ v ] = multiplicative_schwarz_method( A,c,I1,I2,m, a,b, v0,tol,nmax )
%MULTIPLICATIVE_SCHWARZ_METHOD Solves Av=c via multiplicative Schwarz
%method
%   Input: Discretization A,c
%          Restriction operators I1:R->R1, I1:R->R2
%          Initial guess v0, tolerance tol, its max nmax
%          system info m, a,b

h = 1/(2^m);
v_old = v0;
for n=0:nmax
   res_old = c-A*v_old;
   if res_old < tol
       v = v_old;
       break
   end
   %Solve on R1
   ic = 1;
   id = 3;
   %permute f
   f = permute_f_R1(I1*res_old, h,a,b,ic,id);
   %do FFT solve
   solved = rectangular_twod_poissons(m,f',ic,id, zeros(1,id/h-1), ...
       zeros(1,id/h-1),zeros(1,ic/h-1),zeros(1,ic/h-1))';
   %repermute solved output
   right_format = repermute_f_R1(solved, h,a,b);
   
   %compute v^n+1/2
   v_half = v_old + I1'*right_format;
   
   %get next residual
   res_new = c-A*v_half;
   if res_new < tol
       v = v_half;
       break
   end
   
   %Solve on R2
   ic = 4;
   id = 1;
   %permute
   f = permute_f_R2(I2*res_new,h,a,ic,id);
   %do fft solve
   solved = rectangular_twod_poissons(m,f', ic, id, zeros(1,id/h-1), ...
       zeros(1,id/h-1), zeros(1,ic/h-1), zeros(1,ic/h-1))';
   %repermute solved output
   right_format = repermute_f_R2(solved, h,a);
   
   %compute v^n+1
   v_new = v_half + I2'*right_format;
    
   %update
   v_old = v_new;
end

if n==nmax
    v = v_old;
end

end

