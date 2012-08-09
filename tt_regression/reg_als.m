function ttf = reg_als(x, y, fun, r, n)
%REG_ALS Summary of this function goes here
%   d - dimension of the regression function
%   N - number of regression points
%   x - d x N array of coordinates of regression points
%   y - 1 x N array of values of the function to approximate
%   n - d x 1 array, stores number of basis functions in each dimension
%   r - d+1 x 1 array, stores TT-ranks of the tensor of coefficients of the
%       regression function
%   fun - pointer to a function that returns value of the j-th 
%       one-dimensional basis function in an i-th dimension at a point, 
%       syntax: fun(i, j, x).
%   basis - d x 1 cell array of values of basis functions at given points 
%       each element is an (n(i) x N) array
%   fixed - d x 1 cell of aggregated fixed dimensions of the regression 
%   function at k-th step of the ALS algorithm
%   coeff - TT-tensor of coefficient of the regression function

d = numel(n);
coeff = tt_rand(n, d, r);
coeff_ps = coeff.ps;
coeff_cr = coeff.core;

% evaluate basis functions at given points
N = size(x,2);
basis_ps = cumsum([1 ; n*N]);
basis_cr = zeros(basis_ps(d+1) - basis_ps(1), 1);
for dim = 1: d
    basis_ofs = basis_ps(dim);
    n_rows = n(dim);
    for pt = 1: N
        for j = 1: n(dim)
            basis_cr(basis_ofs + (pt-1)*n_rows + (j-1)) = fun(dim, j, x(dim,pt));
        end
    end
end
     
% assemble fixed part of the TT-function
fixed_ps = cumsum([1 ; N ; r(2:d)*N ; N]);
fixed_cr = zeros(fixed_ps(d+2) - fixed_ps(1), 1);
fixed_cr(fixed_ps(d+1): fixed_ps(d+2)-1) = ones(N, r(d+1));
for dim = d: -1: 2
    c = reshape(coeff_cr(coeff_ps(dim): coeff_ps(dim+1)-1), [r(dim) n(dim) r(dim+1)]);
    c = permute(c, [1 3 2]);
    b = reshape(basis_cr(basis_ps(dim): basis_ps(dim+1)-1), [n(dim) N]);
    t1 = ten_conv(c, 3, b); % r(dim) x r(dim+1) x N
    t2 = reshape(fixed_cr(fixed_ps(dim+1): fixed_ps(dim+2)-1), [r(dim+1), N]);
    t3 = zeros(r(dim), N);
    for pt = 1: N
        t3(:,pt) = t1(:,:,pt)*t2(:,pt);
    end
    fixed_cr(fixed_ps(dim): fixed_ps(dim+1)-1) = reshape(t3, [r(dim)*N, 1]);
end
fixed_cr(fixed_ps(1): fixed_ps(2)-1) = ones(N, r(1));

% initial sum of squares
VAR = var(y);

n_swp = 0;
R2 = 0;
swp_tp = 'f';
dim = 1;
while R2 < 0.99999 && n_swp < 10
    n_i = r(dim);
    n_j = n(dim);
    n_k = r(dim+1);
    
    u_i = reshape(fixed_cr(fixed_ps(dim): fixed_ps(dim+1)-1), [N n_i]);
    u_j = reshape(basis_cr(basis_ps(dim): basis_ps(dim+1)-1), [n_j N]);
    u_k = reshape(fixed_cr(fixed_ps(dim+1): fixed_ps(dim+2)-1), [n_k N]);
    
    A = zeros(n_i, n_j, n_k, n_i, n_j, n_k);
    b = zeros(n_i, n_j, n_k);
    
    tic;
    for pt = 1: N
        for i = 1: n_i
            for j = 1: n_j
                for k = 1: n_k
                    b(i,j,k) = b(i,j,k) + u_i(pt,i) * u_j(j,pt) * u_k(k,pt) * y(pt);
                    for i_ = 1: n_i
                        for j_ = 1: n_j
                            for k_ = 1: n_k
                                A(i,j,k,i_,j_,k_) = A(i,j,k,i_,j_,k_) + u_i(pt,i)*u_i(pt,i_) * u_j(j,pt)*u_j(j_,pt) * u_k(k,pt)*u_k(k_,pt);
                            end
                        end
                    end
                end
            end
        end
    end
    toc
    
    A = reshape(A, [n_i*n_j*n_k n_i*n_j*n_k]);
    b = reshape(b, [n_i*n_j*n_k, 1]);
    coeff_cr(coeff_ps(dim):coeff_ps(dim+1)-1) = A\b;
    
    if (swp_tp == 'f' && dim < d) || (swp_tp == 'b' && dim == 1)
        c = reshape(coeff_cr(coeff_ps(dim):coeff_ps(dim+1)-1), [n_i n_j n_k]);
        c = permute(c, [1 3 2]);
        t1 = ten_conv(c, 3, u_j);
        t2 = zeros(N, n_k);
        for pt = 1: N
            t2(pt,:) = u_i(pt,:)*t1(:,:,pt);
        end
        fixed_cr(fixed_ps(dim+1): fixed_ps(dim+2)-1) = reshape(t2, [N*n_k, 1]);
        
        %R^2
        if swp_tp == 'b'
            MSE = 0;
            for pt = 1: N
            MSE = MSE + (t2(pt,:)*u_k(:,pt) - y(pt))^2;
            end
            MSE = MSE/N;
            R2 = 1 - MSE/VAR;
        end
    elseif (swp_tp == 'b' && dim > 1) || (swp_tp == 'f' && dim == d)
        c = reshape(coeff_cr(coeff_ps(dim):coeff_ps(dim+1)-1), [n_i n_j n_k]);
        c = permute(c, [1 3 2]);
        t1 = ten_conv(c, 3, u_j);
        t2 = zeros(n_i, N);
        for pt = 1: N
            t2(:,pt) = t1(:,:,pt)*u_k(:,pt);
        end
        fixed_cr(fixed_ps(dim): fixed_ps(dim+1)-1) = reshape(t2, [n_i*N, 1]);
        
        %R^2
        if swp_tp == 'f'
            MSE = 0;
            for pt = 1: N
                MSE = MSE + (u_i(pt,:)*t2(:,pt) - y(pt))^2;
            end
            MSE = MSE/N;
            R2 = 1 - MSE/VAR;
        end
    end
    
    if swp_tp == 'f'
        if dim < d
            dim = dim + 1;
        elseif dim == d
            swp_tp = 'b';
            dim = dim - 1;
            n_swp = n_swp + 1;
            fprintf('=tt_reg_als= Sweep %d, R^2: %.2f%%, MSE: %1.2e\n', n_swp, 100*R2, MSE);
        end
    elseif swp_tp == 'b'
        if dim > 1
            dim = dim - 1;
        elseif dim == 1
            swp_tp = 'f';
            dim = dim + 1;
            n_swp = n_swp + 1;
            fprintf('=tt_reg_als= Sweep %d, R^2: %.2f%%, MSE: %1.2e\n', n_swp, 100*R2, MSE);
        end
    end
end

coeff.core = coeff_cr;
ttf = tt_function(fun, coeff);

end
