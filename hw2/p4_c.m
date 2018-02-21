%Problem 4c- Hw 2
%Carter Johnson

%Solve LUx = b via Lc = b, Ux = c

%for small system
load('small_ex.mat');
[JL, IL, VL] = coo_to_csc_strictly_lower_tri(L);
c = unit_lower_tri_system_solve(JL, IL, VL, b);
[JU, IU, VU] = coo_to_csc_upper_tri(U);
x = upper_tri_system_solve(JU, IU, VU, c);
x

%Solve large system
tic
load('large_ex.mat');
[JL, IL, VL] = coo_to_csc_strictly_lower_tri(L);
c = unit_lower_tri_system_solve(JL, IL, VL, b);
[JU, IU, VU] = coo_to_csc_upper_tri(U);
x = upper_tri_system_solve(JU, IU, VU, c);
%print out entries
x(50000), x(100000), x(150000), x(200000), x(250000)
toc