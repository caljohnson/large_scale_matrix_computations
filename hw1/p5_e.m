%HW 1 Problem 5e  - MAT226B WQ1028
%Carter Johnson

%given Toeplitz matrix
T = [ 2 -1 4 -3;
      1 2 -1 4;
      -5 1 2 -1;
      6 -5 1 2;];
   
%given vector x
x = [-1; 2; 1; 4;];
y1 = T*x

%construct row vector t to represent T
t = [6 -5 1 2 -1 4 -3];
y2 = toeplitz_matrix_product(t,x)

y1-y2

%now test on n=500,000 case
load('t_and_x.mat')
y = toeplitz_matrix_product(t,x);

%display
display('n=500000 case')
y(1)
y(100000)
y(200000)
y(300000)
y(400000)
y(500000)
sum(y)
norm(y,2)
