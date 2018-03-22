function [ f ] = repermute_f_R2( f_half,h, a )
%REPERMUTE_F_R2 takes f in rectangular shape from FFT solver
% and puts it into vector form f = f_r1\r2, f_gamma2, f_r1_int_r2

%intersection
f_intersect = reshape(f_half(1:(1-a)/h-1, 1:end),[],1);

%gamma boundary
f_gamma = reshape(f_half((1-a)/h, 1:end), [], 1);

%interior box
f_interior = reshape(f_half((1-a)/h+1:end, 1:end), [], 1);

f = [f_intersect; f_gamma; f_interior];

end

