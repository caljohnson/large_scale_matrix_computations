function [ v, H ] = arnoldi_process( A, r, kmax )
%ARNOLDI_PROCESS Constructs the orthonormal basis vectors v1, v2, 
%                                               ..., vkmax+1 for K(k,A,r)
%  INPUT: r - col. vector in C^n
%         A - routine to compute q=Ap
%         kmax - max iterations
%
%  OUTPUT: orthonormal basis v_1, ..., v_kmax+1
%          Upper-hessenberg matrix H_k in C^(kmax+1 x kmax)

v = zeros(size(r,1), kmax+1);
H = zeros(kmax+1,kmax);
beta = norm(r,2);
v(:,1) = r./beta;

for k=1:kmax
   q = A(v(:,k));
   for j=1:k
      H(j,k)  = v(:,j)'*q;
      q = q - H(j,k)*v(:,j);
   end
   H(k+1,k) = norm(q,2);
   if H(k+1,k)==0
       return
   else
       v(:,k+1) = q./H(k+1,k);
   end
end
   
end

