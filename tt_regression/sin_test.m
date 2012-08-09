function ttf = sin_test()

d = 10;
r = [1; 2*ones(d-1,1); 1];
n = 2*ones(d,1);
N = 20000;
x = 2*pi*rand(d,N);
y = zeros(1,N);
for i = 1: N
    y(i) = sin(sum(x(:,i)));
end

ttf = reg_als(x, y, @fun, r, n)

end

function fx = fun(i, j, x)
% i-th dimension, j-th basis function at x

if j == 1
    fx = sin(x);
elseif j == 2
    fx = cos(x);
end

end
