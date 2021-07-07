% Solve heat equation with QTT


% Set d0t=0 for simple time-stepping solution

d0t = 12; % quantics dims for t
d0x = 8; % quantics dims for x
dpx = 2; % phys. dims for x

a = 0;
b = 1;
h = (b-a)/(2^d0x+1);

tol = 1e-5;
eps = 1e-8;

tranges = 0:0.2:0.6;

% negative Laplacian in space
Ax = tt_qlaplace_dd(d0x*ones(1,dpx));
Ax = Ax/(h^2);
Ix = tt_eye(2, dpx*d0x);

% Initial state
x = (a+h:h:b-h)';
u00 = exp(-(x-(a+b)*0.5).^2*16/0.5);
u00 = tt_tensor(reshape(u00, 2*ones(1,d0x)), eps);

u0 = [];
for i=1:dpx
    u0 = tkron(u00,u0);
end

% Fixed source term
f = tt_ones(2, d0x*dpx);

u = u0;
ttimes = zeros(1, max(size(tranges))-1);  % CPU times
resids = zeros(1, max(size(tranges))-1);  % full residuals (space+time)
sp_resids = zeros(1, max(size(tranges))-1); % stabilisation residual (Ax-f)
for t=1:max(size(tranges))-1
    tau = (tranges(t+1)-tranges(t))/(2^d0t);
    
    % Crank-Nicolson matrices
    if (d0t>0)
        Grad_t = IpaS(d0t, -1)/tau;
        CN_t = IpaS(d0t, 1)*0.5;
        e1 = tt_unit(2, d0t, 1);
        et = tt_ones(2, d0t);
    else
        % Just a single time step, but mode size should be >1
        Grad_t = tt_eye(2,1)/tau;
        CN_t = tt_eye(2,1)*0.5;
        e1 = tt_ones(2,1);
        et = tt_ones(2,1);
    end
    
    M = tkron(Ix, Grad_t) + tkron(Ax, CN_t);
    
    tic;
    % Right hand side of global space-time system
    u_rhs = u/tau - (Ax*u)*0.5; % stuff u0 to rhs of CN scheme
    u_rhs = round(u_rhs, eps);
    rhs = tkron(u_rhs, e1) + tkron(f, et);
    rhs = round(rhs, eps);
    % Expand initial state in time
    U = tkron(u, et);
   
    % Solve the global system
    U = amen_solve2(M, rhs, tol, 'x0', U);
    ttimes(t) = toc;    
    resids(t) = norm(M*U-rhs)/norm(rhs);
    
    % extract the last snapshot to prepare a new start
    if (d0t>0)
        ext = tt_unit(2,d0t,2*ones(d0t,1));
    else
        ext = tt_unit(2,1,2);
    end
    u = dot(ext, U, dpx*d0x+1, U.d);
    
    sp_resids(t) = norm(Ax*u-f)/norm(f);
    
    % Plot a 2D marginal
    u2 = tt_reshape(u, 2^d0x*ones(1,dpx));
    if (dpx>2)
        ext = tt_ones(2^d0x, dpx-2);
        u2 = dot(ext, u2, 3, dpx);
    end
    
    mesh(full(u2, 2^d0x*[1,1]));
    drawnow;
    
    fprintf('\nTime range: [%g; %g], solve_time = %g, res = %3.1e, spat_res=%3.3e\n', tranges(t), tranges(t+1), sum(ttimes(1:t)), resids(t), sp_resids(t));    
end
