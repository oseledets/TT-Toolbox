n = 1000; 
m = 1000;

f = @(i,j) 1.0/(i + j);

[u, v] = cross2d(f, n, m, 1e-6);

n = 100;
m = 100;


a = zeros(n, m);

for i = 1:n
    for j = 1:m
        a(i, j) = 1.0/(i + j);
    end
end

[u, v] = cross2d_mat(a, 1e-6);