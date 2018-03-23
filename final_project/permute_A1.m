function [ A1 ] = permute_A1( A1, h, a, b )
%PERMUTE_A1 takes rectangular form A1
%   and puts it in rearranged form to match 
%       with f= f_r1\r2, f_gamma2, f_r1_int_r2


no_points = size(A1,1);

%get new order numbers
if b~=0
    lower_box_nos = 1:(b/h-1)*(1/h-1);
    side_box_nos = lower_box_nos(end)+1:lower_box_nos(end)+(1/h+1)*(a/h-1);
else
    side_box_nos = 1:(1/h+1)*(a/h-1);
end
if b~=2
    top_box_nos = side_box_nos(end)+1:side_box_nos(end)+((3-b-1)/h-1)*(1/h-1);
    lower_gamma_nos = top_box_nos(end)+1:top_box_nos(end)+(1-a)/h;
else
    lower_gamma_nos = side_box_nos(end)+1:side_box_nos(end)+(1-a)/h;
end
side_gamma_nos = lower_gamma_nos(end)+1:lower_gamma_nos(end) + (1/h-1);
if b~=2
    top_gamma_nos = side_gamma_nos(end)+1:side_gamma_nos(end) + (1-a)/h;
    intersect_nos = top_gamma_nos(end)+1:no_points;
else
    intersect_nos = lower_gamma_nos(end)+1:no_points;
end
    
%put them in old order
width1 = a/h-1;
if b~=0
    order = [lower_box_nos side_box_nos(1:width1) lower_gamma_nos];
    if b~=2
        width2 = (1-a)/h-1;
    else
        width2 = (1-a)/h-2;
    end
else
    order = [side_box_nos(1:width1) lower_gamma_nos];
    width2 = (1-a)/h-2;
end

for ii=1:1/h-1
   j = width1*ii+1;
   k = width2*(ii-1)+1;
   order = [order side_box_nos(j:j+width1-1)];
   order = [order side_gamma_nos(ii)];
   order = [order intersect_nos(k:k+width2-1)];
end
if b~=2
    order = [order side_box_nos(end-width1+1:end) top_gamma_nos top_box_nos];
else
    order = [order side_box_nos(end-width1+1:end) side_gamma_nos(end) intersect_nos(end-width2:end)];
end

A1 = A1(order, order);

end

