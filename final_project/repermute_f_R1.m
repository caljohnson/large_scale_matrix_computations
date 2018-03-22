function [ f ] = repermute_f_R1( f_half,h, a, b )
%REPERMUTE_F_R1 takes f in rectangular shape from FFT solver
% and puts it into vector form f = f_r1\r2, f_gamma2, f_r1_int_r2

%bottom interior box
f_int_low = reshape(f_half(1:end, 1:b/h-1), [], 1);

%left interior box
f_int_side = reshape(f_half(1:a/h-1, b/h:(b+1)/h), [], 1);

%top interior box
f_int_hi = reshape(f_half(1:end, (b+1)/h+1:end),[],1);

%bottom gamma2
f_gamma_lo = reshape(f_half(a/h:end,b/h), [],1);

%left gamma2
f_gamma_side = reshape(f_half(a/h,b/h+1:(b+1)/h-1), [], 1);

%top gamma 2
f_gamma_top = reshape(f_half(a/h:end,(b+1)/h), [], 1);

%intersection
f_intersect = reshape(f_half(a/h+1:end, b/h+1:(b+1)/h-1), [], 1);

f = [f_int_low; f_int_side; f_int_hi; f_gamma_lo; ...
    f_gamma_side; f_gamma_top; f_intersect;];

end

