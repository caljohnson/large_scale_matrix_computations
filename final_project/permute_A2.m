function [ A2 ] = permute_A2( A2, h, a )
%PERMUTE_A2 takes rectangular form A2
%   and puts it in rearranged form to match 
%       with f= f_r1_int_r2, f_gamma1, f_r2\r1


no_points = size(A2,1);

%get new order numbers
intersect_nos = 1:((1-a)/h-1)*(1/h-1);
gamma_nos = intersect_nos(end)+1:intersect_nos(end)+(1/h-1);
box_nos = gamma_nos(end)+1:no_points;


%put them in old order
order = [];
width1 = (1-a)/h-1;
width2 = 3/h-1;
for ii=1:1/h-1
   j = width1*(ii-1)+1;
   k = width2*(ii-1)+1;
   order = [order intersect_nos(j:j+width1-1) gamma_nos(ii) box_nos(k:k+width2-1)];
end

A2 = A2(order, order);
end

