function [ w ] = fft_matrix_vector_prod( v,h )
%FFT_MATRIX_VECTOR_PROD Computes w = Zv
% via fft
%   INPUT: vector v
%          grid-spacing h

%obtain vector size m
m = size(v,1);

%pad the input vector
v = [0; v; zeros(m+1,1);];

%apply DFT matrix
w = fft(v);

%obtain matrix-vector product w
w = -sqrt(2*h)*imag(w(2:m+1));

end

