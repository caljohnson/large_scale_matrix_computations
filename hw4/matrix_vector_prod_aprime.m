function [ q ] = matrix_vector_prod_aprime( v, D, D_0, F, G)
%MATRIX_VECTOR_PROD_APRIME Computes q = A'v
%   using A' =  D( (D-G)^{-1} + (D-F)^{-1} (I + D_1 (D-G)^{-1})),
% INPUT: v - vector to multiply by 
%        D - sparse diagonal matrix to approximate D_0 in COO format
%        sparse matrices D_0, F, G s.t. A = D_0 - F - G in COO format

%Compute b = (D-G)^{-1}v' via a triangular solve with D-G.
[JU, IU, VU] = coo_to_csc_upper_tri(D-G); %convert sparse matrices to CSC
b = upper_tri_system_solve(JU, IU, VU, v);

%Compute c = (I + D_1 (D-G)^{-1})v' = v' + D_1 b, 
%with one multiplication with the diagonal entries of D_1 and one SAXPY.
c = v + (D_0-2*D)*b;

%Compute d = (D-F)^{-1}c via a triangular solve with D-F.
[JL, IL, VL] = coo_to_csc_lower_tri(D-F); %convert sparse matrices to CSC
d = lower_tri_system_solve(JL, IL, VL, c);

%Compute q' = D(b+d) with one SAXPY and 
%one multiplication with the diagonal entries of D.
q = D*(b+d);

end

