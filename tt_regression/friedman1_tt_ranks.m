%Friedman 1
d = 5; a = 0; b = 1; n = 20; N = n^d; h = (b-a)/(n-1);
e = ones(1,n);
x1 = a:+h:b;
x = zeros(d,N);
x(1,:) = kron(kron(kron(kron(e,e),e),e),x1);
x(2,:) = kron(kron(kron(kron(e,e),e),x1),e);
x(3,:) = kron(kron(kron(kron(e,e),x1),e),e);
x(4,:) = kron(kron(kron(kron(e,x1),e),e),e);
x(5,:) = kron(kron(kron(kron(x1,e),e),e),e);
y = 10*sin(pi*x(1,:).*x(2,:)) + 20*(x(3,:) - 0.5*ones(1,N)).^2 + 10*x(4,:) + 5*x(5,:);
y = reshape(y, n*ones(1,d));
tt = tt_tensor(y, 1.0e-3)

%Based on observation of the true tt_ranks we set up r and n parameters as follws:
tt_rank = [1 ; 4 ; 2 ; 2 ; 2 ; 1];
n_basis = [5 ; 5 ; 3 ; 2 ; 2];

N = 20000;
x = rand(d,N);
y = 10*sin(pi*x(1,:).*x(2,:)) + 20*(x(3,:) - 0.5*ones(1,N)).^2 + 10*x(4,:) + 5*x(5,:);

fun = @(i, j, x) x^(j-1);
basis_ps = cumsum([1 ; n_basis*N]);
basis_cr = zeros(basis_ps(d+1) - basis_ps(1), 1);
for dim = 1: d
    basis_ofs = basis_ps(dim);
    n_rows = n_basis(dim);
    for pt = 1: N
        for j = 1: n_basis(dim)
            basis_cr(basis_ofs + (pt-1)*n_rows + (j-1)) = fun(dim, j, x(dim,pt));
        end
    end
end

coeff = reg_als(basis_cr, y, tt_rank, n_basis);
ttf = tt_function(fun, coeff);