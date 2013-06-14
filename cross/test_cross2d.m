n = 500;
m = 500;
f = @(i,j) 1.0/(i + 2 * j);
%f = @(i,j) 1.0;
[u, v] = cross2d_new(f, n, m, 1e-12);
