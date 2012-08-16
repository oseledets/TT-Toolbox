function ttf = sin_test()

d = 10;
r = 2;
n = 2;
N = 20000;
x = 2*pi*rand(d,N);
y = sin(sum(x));
tt_rank = [1 ; r*ones(d-1,1) ; 1];
n_basis = n*ones(d,1);
basis_ps = cumsum([1 ; n_basis*N]);
basis_cr = zeros(basis_ps(d+1) - basis_ps(1), 1);
for dim = 1: d
    t = zeros(n, N);
    t(1,:) = sin(x(dim,:));
    t(2,:) = cos(x(dim,:));
	basis_cr(basis_ps(dim):basis_ps(dim+1)-1) = reshape(t, [n*N 1]);
end    
coeff = reg_als(basis_cr, y, tt_rank, n_basis);
ttf = tt_function(@fun, coeff);

end

function fx = fun(~, j, x)
% i-th dimension, j-th basis function at x

if j == 1
    fx = sin(x);
elseif j == 2
    fx = cos(x);
end

end
