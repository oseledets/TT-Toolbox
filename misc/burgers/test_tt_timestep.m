dx = 10;
dt = 14;
alpha = 0.1;

a = -10; b = 10;
T = 10;

h = (b-a)/(2^dx-1);
tau = T/(2^dt);

eps = 1e-6;

u0 = exp(-((a:h:b).^2)');
u0 = tt_tensor(full_to_qtt(u0, eps));
f = -exp(-(((a:h:b)+5).^2)')+0.5;
% f = tt_tensor(tt_ones(dx,2))*0;
f = tt_tensor(full_to_qtt(f, eps));


% Neumann lapl
lp = nd_to_full(tt_qlaplace_dd(dx));
% lp(1,1)=1;
lp(2^dx,2^dx)=1;
lp = tt_matrix(full_to_nd(lp, 1e-10, dx));
lp = lp*alpha/(h^2);

% Gradient
A = spdiags(ones(2^dx,1)*[-1,0,1], [-1 0 1], 2^dx, 2^dx);
% A(1,1)=-1;
A(2^dx,2^dx)=1;
A=A/h;
A = tt_matrix(full_to_nd(full(A), 1e-10, dx));

% Mass
M = spdiags(ones(2^dx,1)*[1/3,1/3,1/3], [-1 0 1], 2^dx, 2^dx);
M = tt_matrix(full_to_nd(full(M), 1e-10, dx));

% Quadratic 3-tensor
Q = matmat2quadr(A,M)*0.25;
Q = round(Q, eps);

% Crank-Nicolson "+" part
B = tt_matrix(tt_eye(2,dx))/tau + lp*0.5;
B = round(B, eps);

% Crank-Nicolson "-" part
CNm = tt_matrix(tt_eye(2,dx))/tau - lp*0.5;
CNm = round(CNm, 1e-12);

% timestepping
ttimes = zeros(2^dt,1);
u = u0;
for t=1:2^dt
    t0 = tic;
    % rhs for implicit ts
    ts_rhs = -convten(convten(Q,u,2), u, 2);
    ts_rhs = round(ts_rhs.tt, eps);
    ts_rhs = ts_rhs + mvk2(CNm, u, eps) + f;
    ts_rhs = round(ts_rhs, eps);
    
    % Matrix for implicit ts
    Mu = mvk2(M, u, eps);
    Au = mvk2(A, u, eps);
    ts_B = (diag(Mu)*A + diag(Au)*M)*0.25 + B;
    ts_B = round(ts_B, eps);
        
    u = dmrg_quadr_solve(ts_B,Q,ts_rhs,u,eps,[], [], 10);
    
    ttimes(t)=toc(t0);
    
    resid = convten(convten(Q,u,2), u, 2);
    resid = round(resid.tt, eps);
    resid = resid + mvk2(ts_B, u, eps) - ts_rhs;
    fprintf('Timestep %d (%3.5f) done. dmrg_resid: %3.3e, time: %3.5f, erank: %g\n', t, tau*t, norm(resid)/norm(ts_rhs), ttimes(t), erank(u));
    plot(full(u, 2^dx), '-');
    str_title = sprintf('Timestep %d (%3.5f)', t, tau*t);
    title(str_title);    
    pause(0.1);    
end;
