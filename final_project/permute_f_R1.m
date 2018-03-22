function [ f ] = permute_f_R1( f_half, h, a, b, ic, id )
%PERMUTE_F_R1 takes vector form f = f_r1\r2, f_gamma2, f_r1_int_r2
%   and puts it into rectangular shape for FFT solver

f = zeros(ic/h-1, id/h-1);

%bottom interior box
no_int_1 = (1/h-1)*(b/h-1);
f(1:end, 1:b/h-1) = reshape(f_half(1:no_int_1), 1/h-1, b/h-1);

%left interior box
no_int_2 = no_int_1+(1/h+1)*(a/h-1);
f(1:a/h-1, b/h:(b+1)/h) = reshape(f_half(no_int_1+1:no_int_2), a/h-1, 1/h+1);

%top interior box
no_int_3 = no_int_2+(1/h-1)*((3-b-1)/h-1);
f(1:end, (b+1)/h+1:end) = reshape(f_half(no_int_2+1:no_int_3), 1/h-1, (3-b-1)/h-1);

%bottom gamma2
no_int_4 = no_int_3 + (1-a)/h;
f(a/h:end,b/h) = f_half(no_int_3+1:no_int_4);

%left gamma2
no_int_5 = no_int_4 + 1/h - 1;
f(a/h,b/h+1:(b+1)/h-1) = f_half(no_int_4+1:no_int_5)';

%top gamma 2
no_int_6 = no_int_5 + (1-a)/h;
f(a/h:end,(b+1)/h) = f_half(no_int_5+1:no_int_6);

%intersection
f(a/h+1:end, b/h+1:(b+1)/h-1) = reshape(f_half(no_int_6+1:end), (1-a)/h-1, 1/h-1);

end

