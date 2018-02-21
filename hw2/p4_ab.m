%Problem 4ab - HW 2
%Carter Johnson

load('small_ex.mat');
%part a
[JL, IL, VL] = coo_to_csc_strictly_lower_tri(L)
%part b
[JU, IU, VU] = coo_to_csc_upper_tri(U)

load('large_ex.mat');
[JL, IL, VL] = coo_to_csc_strictly_lower_tri(L);
[JU, IU, VU] = coo_to_csc_upper_tri(U);
%part a
I = IL;
I(50000), I(100000), I(150000), I(200000), I(250000)
%part b
I = IU;
I(50000), I(100000), I(150000), I(200000), I(250000)


