% Set d0t=0 for simple time-stepping solution

d0t = 12; % quantics dims for t
d0x = 8; % quantics dims for x
dpx = 2; % phys. dims for x

a = 0;
b = 1;
h = (b-a)/(2^d0x+1);

tol = 1e-5;
eps = 1e-8;

tranges = [0:0.5:0.5];

Ax = tt_qlaplace_dd(d0x*ones(1,dpx));
Ax = Ax/(h^2);
Ix = tt_eye(2, dpx*d0x);

x = (a+h:h:b-h)';
u00 = exp(-(x-(a+b)*0.5).^2*16/0.5);
u00 = tt_tensor(reshape(u00, 2*ones(1,d0x)), eps);

u0 = [];
for i=1:dpx
    u0 = kron(u00,u0);
end;

f = tt_ones(2, d0x*dpx); %*0;

u = u0;
ttimes = zeros(1, max(size(tranges))-1);
resids = zeros(1, max(size(tranges))-1);
sp_resids = zeros(1, max(size(tranges))-1);
for t=1:max(size(tranges))-1
    tau = (tranges(t+1)-tranges(t))/(2^d0t);
    
    if (d0t>0)
        Grad_t = IpaS(d0t, -1)/tau;
        CN_t = IpaS(d0t, 1)*0.5;
        e1 = tt_unit(2, d0t, 1);
        et = tt_ones(2, d0t);
    else
        Grad_t = tt_matrix(1)/tau;
        CN_t = tt_matrix(1)*0.5;
        e1 = tt_tensor(1);
        et = tt_tensor(1);
    end;
    
    M = kron(Ix, Grad_t) + kron(Ax, CN_t);
    
    tic;
    u_rhs = u/tau - (Ax*u)*0.5; % stuff u0 to rhs of KN scheme
    u_rhs = round(u_rhs, eps);
    rhs = kron(u_rhs, e1) + kron(f, et);
    rhs = round(rhs, eps);
    U = kron(u, et);
   
    U = als_adapt_solve(M, rhs, tol, 'x0', U, 'nswp', 40, 'max_full_size', 200, 'kickrank', 3, ...
        'pcatype', 'svd', 'tol2', tol/2, 'ismex', true);
    ttimes(t) = toc;    
    resids(t) = norm(M*U-rhs)/norm(rhs);
    
    % prepare new start
    ind = num2cell([1;2]*ones(1,d0x*dpx+d0t), 1);
    for i=1:d0t
        ind{i+dpx*d0x}=2;
    end;
    if (d0t==0)
        ind{dpx*d0x+1}=1;
    end;
    u = U(ind);
    u = tt_reshape(u, 2*ones(1,dpx*d0x));
    
    sp_resids(t) = norm(Ax*u-f)/norm(f);
    
    u2 = tt_reshape(u, 2^d0x*ones(1,dpx));
    ind = num2cell((1:2^d0x)'*ones(1,dpx), 1);
    for i=3:dpx
        ind{i} = 1;
    end;
    u2 = u2(ind);
    
    mesh(full(u2, 2^d0x*[1,1]));
    drawnow;
    
    fprintf('\nTime range: [%g; %g], solve_time = %g, res = %3.1e, spat_res=%3.3e\n', tranges(t), tranges(t+1), sum(ttimes(1:t)), resids(t), sp_resids(t));    
end;
