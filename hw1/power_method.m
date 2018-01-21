function [ lambda, x ] = power_method( matrix_prod, x, eps, k_max )
%POWER_METHOD Computes the dominant eigenvalue+vector
%   Uses matrix-vector products iteratively to compute the dominant
%   eigenvalue lambda and corr. eigenvector x
%Input:  matrix_prod - function to compute matrix vector products Mx with M
%        x - initial vector x in R^n, x != 0
%        eps - parameter eps>0 for convergence check
%        k_max - iteration limit to safeguard against nonconvergence
%Output: Approximation lambda of the dominant eigenvalue of M
%           and corr. approx eigenvector x

%set lambda to max-magnitude x entry
lambda = x(abs(x)==max(abs(x)));
lambda_old = lambda(1); %ensures we chose one entry

for k=1:k_max
   
   %take matrix-vector product
   x = matrix_prod(x); 
   
   %compute new lambda, set to max-magnitude x entry
   lambda = x(abs(x)==max(abs(x))); 
   lambda_new = lambda(1); %ensures we chose one entry
   
   %set new x
   x = (1/lambda_new)*x;
   
   %check whether successive lambdas within relative tolerance
   if abs(lambda_new - lambda_old)<=eps*abs(lambda_old)
      break 
   end
   
   %check whether eigenvalue residual is small - added to speed up in
   %A_alpha case
   if max(abs(matrix_prod(x) - lambda*x))/max(abs(x)) <= 10^(-15)
       break
   end
end

display('Iteration count:')
k
end

