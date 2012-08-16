function coeff = reg_als(basis_cr, y, r, n)
%REG_ALS Summary of this function goes here
%   d - dimension of the regression function
%   N - number of regression points
%   y - 1 x N array of values of the function to approximate
%   n - d x 1 array, stores number of basis functions in each dimension
%   r - d+1 x 1 array, stores TT-ranks of the tensor of coefficients of the
%       regression function
%   basis_cr - array that stores values of basis functions at regression
%   points
%   fixed - d x 1 cell of aggregated fixed dimensions of the regression 
%   function at k-th step of the ALS algorithm
%   coeff - TT-tensor of coefficient of the regression function

N = size(y, 2);
d = numel(n);
coeff = tt_rand(n, d, r);
coeff_ps = coeff.ps;
coeff_cr = coeff.core;
basis_ps = cumsum([1 ; n*N]);
     
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
    t1 = reshape(t1, [r(dim), r(dim+1)*N]);
    t2 = kron(ones(r(dim),1), reshape(t2, [1, r(dim+1)*N]));
    t3 = t1.*t2;
    t3 = reshape(t3, [r(dim) r(dim+1) N]);
	t3 = sum(t3, 2);
    fixed_cr(fixed_ps(dim): fixed_ps(dim+1)-1) = reshape(t3, [r(dim)*N, 1]);
end
fixed_cr(fixed_ps(1): fixed_ps(2)-1) = ones(N, r(1));

% initial sum of squares
VAR = var(y);

n_swp = 0;
R2 = 0;
swp_tp = 'f';
dim = 1;
while R2 < 0.9999 && n_swp < 10
    n1 = r(dim);
    n2 = n(dim);
    n3 = r(dim+1);
    
    u1 = reshape(fixed_cr(fixed_ps(dim): fixed_ps(dim+1)-1), [N n1]);
    u2 = reshape(basis_cr(basis_ps(dim): basis_ps(dim+1)-1), [n2 N]);
    u3 = reshape(fixed_cr(fixed_ps(dim+1): fixed_ps(dim+2)-1), [n3 N]);
    
    t1 = kron(kron(ones(n3,1), ones(n2,1)), u1');
    t2 = kron(kron(ones(n3,1), u2), ones(n1,1));
    t3 = kron(kron(u3, ones(n2,1)), ones(n1,1));
    t = t1.*t2.*t3;
    
    A = t*t';
    b = t*y';
    c = A\b;
    coeff_cr(coeff_ps(dim):coeff_ps(dim+1)-1) = c;
    
    if (swp_tp == 'f' && dim < d) || (swp_tp == 'b' && dim == 1)
        c = permute(reshape(c, [n1 n2 n3]), [1 3 2]);
        t1 = ten_conv(c, 3, u2);
        t1 = reshape(t1, [n1*n3, N]);
        u1 = kron(ones(1,n3), u1);
        t2 = u1.*t1';
        t2 = reshape(t2, [N n1 n3]);
        t2 = sum(t2, 2);
        t2 = reshape(t2, [N, n3]);
        fixed_cr(fixed_ps(dim+1): fixed_ps(dim+2)-1) = reshape(t2, [N*n3, 1]);
        
        %R^2
        if swp_tp == 'b'
            MSE = sum((sum(t2'.*u3,1) - y).^2)/N;
            R2 = 1 - MSE/VAR;
        end
    elseif (swp_tp == 'b' && dim > 1) || (swp_tp == 'f' && dim == d)
        c = permute(reshape(c, [n1 n2 n3]), [1 3 2]);
        t1 = ten_conv(c, 3, u2);
        t1 = reshape(t1, [n1, n3*N]);
        u3 = kron(ones(n1,1), reshape(u3, [1, n3*N]));
        t2 = t1.*u3;
        t2 = reshape(t2, [n1 n3 N]);
        t2 = sum(t2, 2);
        t2 = reshape(t2, [n1 N]);
        fixed_cr(fixed_ps(dim): fixed_ps(dim+1)-1) = reshape(t2, [n1*N, 1]);
        
        %R^2
        if swp_tp == 'f'
            MSE = sum((sum(u1'.*t2,1) - y).^2)/N;
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

end