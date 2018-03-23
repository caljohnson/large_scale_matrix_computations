function [ f ] = repermute_f_R1( f_half,h, a, b )
%REPERMUTE_F_R1 takes f in rectangular shape from FFT solver
% and puts it into vector form f = f_r1\r2, f_gamma2, f_r1_int_r2

%bottom interior box
if b~=0
    f_int_low = reshape(f_half(1:end, 1:b/h-1), [], 1);
end

%left interior box
if b~=0
    if b~=2
        f_int_side = reshape(f_half(1:a/h-1, b/h:(b+1)/h), [], 1);
    else
        f_int_side = reshape(f_half(1:a/h-1, b/h:3/h-1), [], 1);
    end
else
    f_int_side = reshape(f_half(1:a/h-1, 1:(b+1)/h), [], 1);
end

%top interior box
if b~=2
    f_int_hi = reshape(f_half(1:end, (b+1)/h+1:end),[],1);
end

%bottom gamma2
if b~=0
    f_gamma_lo = reshape(f_half(a/h:end,b/h), [],1);
end

%left gamma2
f_gamma_side = reshape(f_half(a/h,b/h+1:(b+1)/h-1), [], 1);

%top gamma 2
if b~=2
    f_gamma_top = reshape(f_half(a/h:end,(b+1)/h), [], 1);
end

%intersection
f_intersect = reshape(f_half(a/h+1:end, b/h+1:(b+1)/h-1), [], 1);

if b~=0
    if b~=2
        f = [f_int_low; f_int_side; f_int_hi; f_gamma_lo; ...
            f_gamma_side; f_gamma_top; f_intersect;];
    else
         f = [f_int_low; f_int_side; f_gamma_lo; ...
            f_gamma_side; f_intersect;];
    end
else
    f = [f_int_side; f_int_hi; f_gamma_side; f_gamma_top; f_intersect;];
end

end

