function [ VL ] = incomplete_chol( A, J, I )
%INCOMPLETE_CHOL Computes the incomplete Cholesky factorization of A
%   for sparse symmetric matrix A and sparsity structure E
%   where E is specified in CSC format by row-indices J and pointers I by
%       E = {(J(i),k) | i = I(k), I(k)+1, ..., I(k+1)-1, k=1,2,...,n}
%   Input: sparse symmetric matrix A
%           row-indices vector J
%           pointer vector I
%   Output: columnwise-ordered vector VL of Chol factor values

n = size(A,1);
VL = zeros(size(J));

%set l_jk = a_jk for all (j,k) in E
for k=1:n
    VL(I(k):I(k+1)-1) = A(J(I(k):I(k+1)-1),k);
end

for k=1:n
   %check if l_kk <= 0
   kk = I(k); %I(k) points to where col k starts, i.e., row k
   if VL(kk)<=0
       disp('l_kk <= 0, Algorithm fails,');
       return
   end
   %set l_kk = sqrt(l_kk)
   VL(kk) = sqrt(VL(kk));
   %find (j,k) in E w/ j>k
   jks = I(k)+1:I(k+1)-1;
   %set l_jk = (1/l_kk) l_jk for all (j,k) in E w/ j>k
   VL(jks) = VL(jks)/VL(kk);
   %for all (j,k) in E w/ j>k do
   for jk = jks
      j = J(jk);
      %find (i,j) in E
      if j <= n+1
          ij = I(j):I(j+1)-1;
      else
          ij = [];
      end
      %find (i,k)
      temp1 = J(ij);
      temp2 = I(k):I(k+1)-1;
      [~, ij2, ik] = intersect(temp1,J(temp2));
      ij2 = ij(ij2);
      ik = temp2(ik);
      %set l_ij = l_ij - l_jk l_ik for all (i,j) in E
      if ~isempty(ik)
        VL(ij2) = VL(ij2) - VL(jk)*VL(ik);
      end
   end
end

end

