function [D,x]=spectral_diff_pc(N)
h = 2*pi/N;
x = h*(1:N)';
column = [0, 0.5*(-1).^(1:N-1).*cot((1:N-1)*h/2)];
D = toeplitz(column, column([1 N:-1:2]));

end