%HW 1 Problem 5b  - MAT226B WQ1028
%Carter Johnson

%given circulant matrix C
C = [2 -1 0 -3; 
    -3 2 -1 0; 
    0 -3 2 -1; 
    -1 0 -3 2;];
c = [2 -3 0 -1;] %first column of C, defines C
%given vector x
x = [-1; 2; 1; 4;];
%dimension of space
n = 4;

%test matrix-vector product method
y1 = C*x;
y2 = circulant_matrix_product(c,x);

format long e
y1 - y2