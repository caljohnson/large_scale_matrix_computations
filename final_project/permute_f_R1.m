function [ f ] = permute_f_R1( f_half, h, a, b, ic, id )
%PERMUTE_F_R1 takes vector form f = f_r1\r2, f_gamma2, f_r1_int_r2
%   and puts it into rectangular shape for FFT solver

f = zeros(ic/h-1, id/h-1);

%bottom interior box
if b~=0
    no_int_1 = (1/h-1)*(b/h-1);
    f(1:end, 1:b/h-1) = reshape(f_half(1:no_int_1), 1/h-1, b/h-1);
end

%left interior box
if b~=0
    if b~=2
        no_int_2 = no_int_1+(1/h+1)*(a/h-1);
        f(1:a/h-1, b/h:(b+1)/h) = reshape(f_half(no_int_1+1:no_int_2), a/h-1, 1/h+1);
    else
        no_int_2 = no_int_1+(1/h)*(a/h-1);
        f(1:a/h-1, b/h:(b+1)/h-1) = reshape(f_half(no_int_1+1:no_int_2), a/h-1, 1/h);   
    end
else
    no_int_2 = (1/h)*(a/h-1);
    f(1:a/h-1, 1:(b+1)/h) = reshape(f_half(1:no_int_2), a/h-1, 1/h);
end

%top interior box
if b~=2
    no_int_3 = no_int_2+(1/h-1)*((3-b-1)/h-1);
    f(1:end, (b+1)/h+1:end) = reshape(f_half(no_int_2+1:no_int_3), 1/h-1, (3-b-1)/h-1);
else
    no_int_3 = no_int_2;
end

%bottom gamma2
if b~=0
    no_int_4 = no_int_3 + (1-a)/h;
    f(a/h:end,b/h) = f_half(no_int_3+1:no_int_4);
else
    no_int_4 = no_int_3;
end

%left gamma2
no_int_5 = no_int_4 + 1/h - 1;
f(a/h,b/h+1:(b+1)/h-1) = f_half(no_int_4+1:no_int_5)';

%top gamma 2
if b~=2
    no_int_6 = no_int_5 + (1-a)/h;
    f(a/h:end,(b+1)/h) = f_half(no_int_5+1:no_int_6);
else
    no_int_6 = no_int_5;
end

%intersection
if b~=2
    f(a/h+1:end, b/h+1:(b+1)/h-1) = reshape(f_half(no_int_6+1:end), (1-a)/h-1, 1/h-1);
else
    f(a/h+1:end, b/h+1:end) = reshape(f_half(no_int_6+1:end), (1-a)/h-1, 1/h-1);
end

end

