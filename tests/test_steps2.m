%Crank-Nicolson scheme with local time-stepping
d0t = 10; % quantics dims for t
d0x = 9; % quantics dims for x
dpx = 2; % phys. dims for x

a=0; 
b=1;
h = (b-a)/(2^d0x+1);

tol = 1e-6;
eps = 1e-8;
maxit=20;

T = 0.5;
tau = T/(2^d0t);

Ax = tt_qlaplace_dd(d0x*ones(1,dpx));
Ax = Ax/(h^2);

x = (-a+h:h:b-h)';
u0 = kron(sin(pi*x), sin(pi*x));
x = tt_tensor(x);
x = tt_reshape(x, 2*ones(d0x, 1), eps);
e1 = tt_ones(2, d0x);

u0 = tt_reshape(tt_tensor(u0), 2*ones(1,2*d0x), eps);
% u0 = kron(x.*(e1-x), e1).*kron(e1, x.*(e1-x));
% u0 = round(u0, eps);
% u0 = tt_tensor(tt_ones(d0x*dpx,2));

spacial_rhs = u0*0;


% Now, generate parts of KN scheme:
KNm = tt_eye(2,dpx*d0x)/tau - Ax*0.5;
KNm = round(KNm, eps);

KNp = tt_eye(2,dpx*d0x)/tau + Ax*0.5;
KNp = round(KNp, eps);

Euler = tt_eye(2,dpx*d0x) - Ax*tau;
Euler = round(Euler, eps);

ttimes = zeros(1,2^d0t);
resids = zeros(1,2^d0t);
dmrg_resids = zeros(1,2^d0t);
dmrg_swps = zeros(1,2^d0t);
eranks = zeros(1,2^d0t);
u = u0;
for t=1:1:2^d0t
    times_0 = tic;    
    u_new = mvk3(KNm, u, tol, 'nswp', 20, 'y0', u); % 1\tI - A/2
    u_new = u_new + spacial_rhs;
    u = mvk3(Euler, u, tol, 'nswp', 20, 'y0', u) + tau*spacial_rhs;    
%     for i=1:maxit
%         [u,swps] = dmrg_solve2(KNp, u_new, u, tol, [], [], 100, []);
        [u,swps] = dmrg_solve2(KNp, u_new, tol, 'x0', u, 'nswp', 1, 'min_drank', 0, 'ddrank', 0, 'min_dpow', 0);
%         dmrg_resid = norm(KNp*u-u_new)/norm(u_new);
%         fprintf('-- dmrg sweep: %d, resid: %3.3e\n', i, dmrg_resid);
%         if (dmrg_resid<tol*5)
%             break;
%         end;
%     end;
    ttimes(t)=toc(times_0);    
    
    mesh(full(u, 2^d0x*[1,1]))
    appr = norm(u-u0*exp(-2*pi^2*t*tau))/norm(u0*exp(-2*pi^2*t*tau))
    
    dmrg_resid = norm(KNp*u-u_new)/norm(u_new);
    dmrg_swps(t)=swps;
    
    dmrg_resids(t) = dmrg_resid;
    resids(t) = norm(Ax*u - spacial_rhs)/norm(spacial_rhs);
%     resids(t) = dot(Ax*u, u)/norm(Ax*u)/norm(u)-1;
    eranks(t) = erank(u);
    
    fprintf('\nTime step %d (%3.5e) done. ttimes: %3.5f, dmrg_resid: %3.3e, sweeps: %d, spac_resid: %3.5e, erank: %3.2f\n', t, t*tau, ttimes(t), dmrg_resids(t), dmrg_swps(t), resids(t), eranks(t));
%     mesh(full(u, [2^d0x,2^d0x]));    
    pause(0.01);
%     keyboard;
end;
