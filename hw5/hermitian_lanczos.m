function [ T ] = hermitian_lanczos( A,r, kmax )
%HERMITIAN_LANCZOS implements the Hermitian Lanczos process for Hermitian A
% with vector r!=0 to create tridiagonal matrix Tk after k steps
%   INPUT: A - matrix %routine to compute q=Av
%          r - col. vector, initial residual, size nx1
%          kmax - max its size
%   OUTPUT: T - tridiagonal matrix of size kmax x kmax
%           

%allocate space
beta = zeros(kmax+2,1);
alpha = zeros(kmax+1,1);

%Start algorithm
beta(1) = norm(r,2);
%only keep track of 2 V vectors, not full matrix
v_new = r./beta(1);

for k=1:kmax+1
    q = A*v_new;
    if k>1
        q = q - beta(k)*v_old;
    end
    alpha(k) = v_new'*q;
    q = q - alpha(k)*v_new;
    beta(k+1) = norm(q,2);
    if beta(k+1)==0 %k = d(A,r), done
%         T = spdiags([beta(2:k+1),alpha(1:k),beta(2:k+1)], -1:1, kmax, kmax);
        display('k=d(A,r)');
        return
    end
    v_old = v_new;
    v_new = q./beta(k+1);
end

T = full(spdiags([beta(2:end),alpha,beta(1:end-1)], -1:1, kmax,kmax));

end

