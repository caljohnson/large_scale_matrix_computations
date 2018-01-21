%Hw 1 Problem 2 - MAT226B WQ1028
%Carter Johnson

%outgoing-links row-stochastic matrix A
A = [ 0 0 0 0 1/2 0 1/2 0 0 0; 
0 1/4 0 0 0 1/4 1/4 0 1/4 0; 
1/10 1/10 1/10 1/10 1/10 1/10 1/10 1/10 1/10 1/10;
0 0 0 1/2 0 0 1/2 0 0 0;
0 0 0 0 0 0 1/3 0 1/3 1/3;
1/10 1/10 1/10 1/10 1/10 1/10 1/10 1/10 1/10 1/10;
1/10 1/10 1/10 1/10 1/10 1/10 1/10 1/10 1/10 1/10;
0 0 1/5 1/5 0 0 1/5 1/5 1/5 0;
0 1/5 0 1/5 0 0 1/5 1/5 0 1/5;
1/6 0 1/6 1/6 0 0 1/6 0 1/6 1/6;];

%pagerankings is the eigenvector of A' with eigvalue 1
[V, D] = eig(A');

D(1) %eigenvalue 1
x = V(:,1); %column is the right-eigenvector corresponding to eigenvalue 1
x %PageRank vector