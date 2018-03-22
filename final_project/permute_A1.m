function [ A1 ] = permute_A1( A1, h, a, b )
%PERMUTE_A1 takes rectangular form A1
%   and puts it in rearranged form to match 
%       with f= f_r1\r2, f_gamma2, f_r1_int_r2


no_points = size(A1,1);

%get new order numbers
lower_box_nos = 1:(b/h-1)*(1/h-1);
side_box_nos = lower_box_nos(end)+1:lower_box_nos(end)+(1/h+1)*(a/h-1);
top_box_nos = side_box_nos(end)+1:side_box_nos(end)+((3-b-1)/h-1)*(1/h-1);
lower_gamma_nos = top_box_nos(end)+1:top_box_nos(end)+(1-a)/h;
side_gamma_nos = lower_gamma_nos(end)+1:lower_gamma_nos(end) + (1/h-1);
top_gamma_nos = side_gamma_nos(end)+1:side_gamma_nos(end) + (1-a)/h;
intersect_nos = top_gamma_nos(end)+1:no_points;

%put them in old order
width1 = a/h-1;
order = [lower_box_nos side_box_nos(1:width1) lower_gamma_nos];
width2 = (1-a)/h-1;
for ii=1:1/h-1
   j = width1*ii+1;
   k = width2*(ii-1)+1;
   order = [order side_box_nos(j:j+width1-1) side_gamma_nos(ii) intersect_nos(k:k+width2-1)];
end
order = [order side_box_nos(end-width1+1:end) top_gamma_nos top_box_nos];

A1 = A1(order, order);
end

