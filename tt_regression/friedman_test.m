function ttf = friedman_test(n_problem, r, n)
%FRIEDMAN regression

fun = @(i, j, x) x^(j-1);

if n_problem == 1
    %Friedman 1
    d = 10;
    N = 20000;
    x = rand(d,N);
    y = 10*sin(pi*x(1,:).*x(2,:)) + 20*(x(3,:) - 0.5*ones(1,N)).^2 + 10*x(4,:) + 5*x(5,:);
    
    tt_rank = [1 ; r*ones(d-1,1) ; 1];
    n_basis = n*ones(d,1);
    basis_ps = cumsum([1 ; n_basis*N]);
    basis_cr = zeros(basis_ps(d+1) - basis_ps(1), 1);
    for dim = 1: d
        t = zeros(n, N);
        for i = 1: n
            t(i,:) = x(dim,:).^(i-1);
        end
        basis_cr(basis_ps(dim):basis_ps(dim+1)-1) = reshape(t, [n*N 1]);
    end
elseif n_problem == 2
    %Friedman 2
    d = 4;
    N = 20000;
    x = rand(d,N);
    x1 = x(1,:)*100;
    x2 = x(2,:)*520*pi + 40*pi*ones(1,N);
    x3 = x(3,:);
    x4 = x(4,:)*10 + ones(1,N);
    y = sqrt(x1.^2 + (x2.*x3 - 1./(x2.*x4)).^2);
    
    tt_rank = [1 ; r*ones(d-1,1) ; 1];
    n_basis = n*ones(d,1);
    basis_ps = cumsum([1 ; n_basis*N]);
    basis_cr = zeros(basis_ps(d+1) - basis_ps(1), 1);
    for dim = 1: d
        t = zeros(n, N);
        for i = 1: n
            t(i,:) = x(dim,:).^(i-1);
        end
        basis_cr(basis_ps(dim):basis_ps(dim+1)-1) = reshape(t, [n*N 1]);
    end
else
    %Friedman 3
    d = 4;
    N = 20000;
    x = rand(d,N);
    x1 = x(1,:)*100;
    x2 = x(2,:)*520*pi + 40*pi*ones(1,N);
    x3 = x(3,:);
    x4 = x(4,:)*10 + ones(1,N);
    y = atan((x2.*x3 - 1./(x2.*x4))./x1);
    
    tt_rank = [1 ; r*ones(d-1,1) ; 1];
    n_basis = n*ones(d,1);
    basis_ps = cumsum([1 ; n_basis*N]);
    basis_cr = zeros(basis_ps(d+1) - basis_ps(1), 1);
    for dim = 1: d
        t = zeros(n, N);
        for i = 1: n
            t(i,:) = x(dim,:).^(i-1);
        end
        basis_cr(basis_ps(dim):basis_ps(dim+1)-1) = reshape(t, [n*N 1]);
    end
end

coeff = reg_als(basis_cr, y, tt_rank, n_basis);
ttf = tt_function(fun, coeff);

end

