function [ y ] = twod_fft_matrix( x )
%TWOD_FFT_MATRIX Computes X' = Z'XZ via fft


function y = oned_fft_matrix(x) 
n = size(x,1); 
m = size(x,2);
y = [zeros(1,m);x];
y = imag(fft(y,2*n+2)); 
y = sqrt(2/(n+1))*y(2:n+1,:);
end

y = oned_fft_matrix(oned_fft_matrix(x)')';

end

