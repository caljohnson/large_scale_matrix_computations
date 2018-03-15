%Problem 1 - HW 4
%MAT 226B


n=5
I = eye(n,n)
a = 2*ones(5,1)
A = [zeros(1,5) 1; I a;]
b = [1;0;0;0;0;0;];

for i =0:n
    A^i*b
end