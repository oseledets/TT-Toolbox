d0t = 9; % quantics dims for t
d0x = 8; % quantics dims for x
dpx = 2; % phys. dims for x

a=0; 
b=1;
h = (b-a)/(2^d0x+1);

tol = 1e-6;
eps = 1e-8;
maxit=20;

T = 1;
tau = T/(2^d0t);

Ax = tt_matrix(tt_qlaplace_dd(d0x*ones(1,dpx)));
Ax = Ax/(h^2);

u0 = tt_tensor(tt_ones(d0x*dpx,2));
spacial_rhs = u0;


% Now, generate parts of KN scheme:
KNm = tt_matrix(tt_eye(2,dpx*d0x))/tau - Ax*0.5;
KNm = round(KNm, eps);

KNp = tt_matrix(tt_eye(2,dpx*d0x))/tau + Ax*0.5;
KNp = round(KNp, eps);

Euler = tt_matrix(tt_eye(2,dpx*d0x)) - Ax*tau;
Euler = round(Euler, eps);

ttimes = zeros(1,2^d0t);
resids = zeros(1,2^d0t);
dmrg_resids = zeros(1,2^d0t);
dmrg_swps = zeros(1,2^d0t);
eranks = zeros(1,2^d0t);
u = u0;
for t=1:1:2^d0t
    times_0 = tic;    
    u_new = mvk(KNm, u, tol, 20, u, 1000); % 1\tI - A/2
    u_new = u_new + spacial_rhs;
    u = mvk(Euler, u, tol, 20, u, 1000) + tau*spacial_rhs;    
%     for i=1:maxit
        [u,swps] = dmrg_solve2(KNp, u_new, u, tol, [], [], 100, []);
%         dmrg_resid = norm(KNp*u-u_new)/norm(u_new);
%         fprintf('-- dmrg sweep: %d, resid: %3.3e\n', i, dmrg_resid);
%         if (dmrg_resid<tol*5)
%             break;
%         end;
%     end;
    ttimes(t)=toc(times_0);    
    
    dmrg_resid = norm(KNp*u-u_new)/norm(u_new);
    dmrg_swps(t)=swps;
    
    dmrg_resids(t) = dmrg_resid;
    resids(t) = norm(Ax*u - spacial_rhs)/norm(spacial_rhs);
    eranks(t) = erank(u);
    
    fprintf('\nTime step %d (%3.5e) done. ttimes: %3.5f, dmrg_resid: %3.3e, sweeps: %d, spac_resid: %3.5e, erank: %3.2f\n', t, t*tau, ttimes(t), dmrg_resids(t), dmrg_swps(t), resids(t), eranks(t));
%     mesh(full(u, [2^d0x,2^d0x]));    
%     pause(0.1);
%     keyboard;
end;
