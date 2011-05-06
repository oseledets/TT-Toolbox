d0t = 20; % quantics dims for t
d0x = 9; % quantics dims for x
dpx = 2; % phys. dims for x

tol = 1e-6;
maxit = 20;

tau = 1e-6;

St = tt_shf(d0t); St = tt_matrix(St); St=St'; % lower shift matrix for gradient
It = tt_matrix(tt_eye(2,d0t));
Grad_t = (It-St)/tau; % for gradient
Grad_t = round(Grad_t, 1e-12);

KN_term = (It+St)*0.5; % Krank-Nikolson term
KN_term = round(KN_term, 1e-12);

e1t = cell(d0t,1);
for i=1:d0t
    e1t{i}=[1;0];
end;
e1t = tt_tensor(e1t); % first identity vector for t - we need to put u0 into rhs


Ax = tt_matrix(tt_qlaplace_dd(d0x*ones(1,dpx)));
Ax = Ax*((2^d0x+1)^2);

% For Fokker-Plank
beta = 1;
ddd = 0.5;
zzz = 0.1;

Sx = tt_shf(d0x);
Grad_x = (Sx - Sx')*(0.5*(2^d0x+1));
Grad_x = round(Grad_x, 1e-10);
Ix = tt_matrix(tt_eye(2, d0x));

Cx1 = kron(Ix, kron(Ix, Grad_x));
Cx2 = kron(Ix, kron(Grad_x, Ix));
Cx3 = kron(Grad_x, kron(Ix, Ix));

ttx = tt_tensor(tt_x(d0x))/(2^d0x+1);
ex = tt_tensor(tt_ones(d0x,2));
x1 = kron(kron(ex,ex), ttx);
x2 = kron(kron(ex,ttx), ex);
x3 = kron(kron(ttx,ex), ex);

r2 = x1.*x1 + x2.*x2 + x3.*x3;
r2 = round(r2, 1e-10);
eexp = funcrs2(r2, @(r)exp(-0.5*r), tol, r2, 20);
v1 = (eexp.*x1)*(z/(2*d^5)) - x1 + beta*x2;
v1 = round(v1, tol);
v2 = (eexp.*x2)*(z/(2*d^5)) - x2;
v3 = (eexp.*x3)*(z/(2*d^5)) - x3;

Ix3 = kron(Ix, kron(Ix, Ix));

% global matrix
M = kron(Ix, Grad_t) + kron(Ax, KN_term);
M = round(M, 1e-12);

% f1 = sin(pi*(1:1:2^d0x)'/(1+2^d0x));
% u0 = full_to_qtt(f1, 1e-12);
% u0 = tt_tensor(u0);
% u0 = kron(u0,u0);

u0 = tt_tensor(tt_ones(d0x*dpx, 2));
u0_rhs = u0/tau - (Ax*u0)*0.5; % stuff u0 to rhs of KN scheme
u0_rhs = round(u0_rhs, 1e-12);

rhs = kron(u0_rhs, e1t);

x = tt_random(2, rhs.d, 2);

results = zeros(maxit,6);
for i=1:maxit
    tic;
    x = dmrg_solve2(M, rhs, x, tol, [], [], 1, M');
    cur_time = toc;
    
    resid_true = norm(M*x-rhs)/norm(rhs);
    resid = norm(M'*(M*x) - M'*rhs)/norm(M'*rhs);
    
    fprintf('\n\n\t cur_time: %g\n\t true_resid: %3.3e\n\t norm.resid: %3.3e\n\t erank: %g\n', cur_time, resid_true, resid, erank(x));
    results(i,1)=i;
    results(i,2)=resid_true;
    results(i,3)=resid;
    results(i,4)=erank(x);
    results(i,5)=cur_time;
    pause(2);
end;

results(:,6)=cumsum(results(:,5));

fprintf('sweep\t true resid\t norm. resid \t erank  \t sw. time \t full time\n');
for i=1:maxit
    fprintf('%d\t %3.3e\t %3.3e\t %3.3f  \t %3.5f \t %3.5f\n', results(i,1),results(i,2),results(i,3),results(i,4),results(i,5),results(i,6));
end;
