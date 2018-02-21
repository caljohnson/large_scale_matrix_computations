%Problem 5 - HW 2
%Carter Johnson

%solve Ax = b using funtions from 4
load('large_ex2.mat'); %load('large_ex1.mat');

%(i)  without column reordering
display('w/o col reorder')
p0 = 1:size(A,1);
%w/o scaling
scaling = 0;
display('w/o scaling:')
x = sparse_system_solve_lu(A,b,p0,scaling);

%w/ scaling
scaling = 1;
display('w/ scaling:')
x = sparse_system_solve_lu(A,b,p0,scaling);

%(ii) w/ column reordering 'colamd(A)'
display('w/ col reorder colamd')
p0 = colamd(A);
%w/o scaling
scaling = 0;
display('w/o scaling:')
x = sparse_system_solve_lu(A,b,p0,scaling);

%w/ scaling
scaling = 1;
display('w/ scaling:')
x = sparse_system_solve_lu(A,b,p0,scaling);

%(iii) w/ column reordering 'colperm(A)'
display('w/ col reorder colperm')
p0 = colperm(A);
%w/o scaling
scaling = 0;
display('w/o scaling:')
x = sparse_system_solve_lu(A,b,p0,scaling);

%w/ scaling
scaling = 1;
display('w/ scaling:')
x = sparse_system_solve_lu(A,b,p0,scaling);



