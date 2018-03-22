function [ f ] = permute_f_R2( f_half, h, a, ic, id )
%PERMUTE_F_R2 takes vector form f = f_r1_int_r2, f_gamma1, f_r2\r1
%   and puts it into rectangular shape for FFT solver

f = zeros(ic/h-1, id/h-1);

%left box, intersection of R1 and R2
no_int_1 = ((1-a)/h-1)*(1/h-1);
f(1:(1-a)/h-1, 1:end) = reshape(f_half(1:no_int_1), (1-a)/h-1, 1/h-1);

%gamma boundary
no_int_2 = no_int_1+1/h-1;
f((1-a)/h, 1:end) = reshape(f_half(no_int_1+1:no_int_2), 1, 1/h-1);

%interior box
f((1-a)/h+1:end, 1:end) = reshape(f_half(no_int_2+1:end), (4-a)/h-1, 1/h-1);
end

