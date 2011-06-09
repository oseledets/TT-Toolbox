d0t = 9; % quantics dims for t
d0x = 9; % quantics dims for x
dpx = 3; % phys. dims for x

a = -10;
b = 10;
% a=0; b=1;
h = (b-a)/(2^d0x+1);

tol = 1e-8;
eps = 1e-8;
maxit = 20;

T = 2;
tau = T/(2^d0t);

Ax = tt_matrix(tt_qlaplace_dd(d0x*ones(1,dpx)));
Ax = Ax/(h^2);

% For Fokker-Plank
beta = 1;
ddd = 0.5;
zzz = 0.1*2;

Sx = tt_matrix(tt_shf(d0x));
Grad_x = (Sx - Sx')/(2*h);
Grad_x = round(Grad_x, eps);
Ix = tt_matrix(tt_eye(2, d0x));

Cx1 = kron(Ix, kron(Ix, Grad_x));
Cx2 = kron(Ix, kron(Grad_x, Ix));
Cx3 = kron(Grad_x, kron(Ix, Ix));

x = a + (0:1:2^d0x-1)'*h;
eexp1 = exp(-0.5*(x.^2)/(ddd^2));
eexp1 = full_to_qtt(eexp1, eps);
eexp1 = tt_tensor(eexp1);

ttx = a*tt_tensor(tt_ones(d0x,2))+tt_tensor(tt_x(d0x,2))*h;
ex = tt_tensor(tt_ones(d0x,2));
x1 = kron(kron(ex,ex), ttx);
x2 = kron(kron(ex,ttx), ex);
x3 = kron(kron(ttx,ex), ex);

eexp = kron(kron(eexp1, eexp1), eexp1);

r2 = x1.*x1 + x2.*x2 + x3.*x3;
% r2 = round(r2, eps);
eexp = funcross(r2, @(r)exp(-0.5*r/(ddd^2)), eps, r2, 10);
% eexp = funcross(x2.*x2, @(r)exp(-0.5*r/(ddd^2)), eps, x2, 10);



% v1 = 0.5*(eexp.*x1)*(zzz/(2*ddd^5)) - x1*0.5;
v1 = 0.5*((eexp.*x1)*(zzz/(2*ddd^5)) - x1);
v1 = round(v1, eps);
v2 = 0.5*((eexp.*x2)*(zzz/(2*ddd^5)) - x2); % - x2.*x2.*x2;
% v3 = 0.5*(eexp.*x3)*(zzz/(2*ddd^5)) - x3*0.5;
v3 = 0.5*((eexp.*x3)*(zzz/(2*ddd^5)) - x3);

vv1 = diag(v1 + beta*x2);
vv2 = diag(v2);
vv3 = diag(v3);

phi = r2*0.5 + 0.5*zzz/(ddd^3)*eexp; % + x2.*x2.*x2.*x2/4;
phi = round(phi, eps);
x_ex = funcross(phi, @(r)exp(-r), eps, phi, 10);

Ax = Ax*0.5 + Cx1*vv1 + Cx2*vv2 + Cx3*vv3;
Ax = round(Ax, eps);

norm_Ax = norm(Ax*x_ex)/norm(x_ex)

u0 = x_ex.*x_ex;
u0 = round(u0, tol);

% f1 = sin((a+(1:1:2^d0x)*h)'*pi);
% f1 = 

% Now, generate parts of KN scheme:
KNm = tt_matrix(tt_eye(2,dpx*d0x))/tau - Ax*0.5;
KNm = round(KNm, eps);

KNp = tt_matrix(tt_eye(2,dpx*d0x))/tau + Ax*0.5;
KNp = round(KNp, eps);

Euler = tt_matrix(tt_eye(2,dpx*d0x)) - Ax*tau;
Euler = round(Euler, eps);

ons = tt_tensor(tt_ones(eexp.d, 2));

results = zeros(2^d0t,1);
angles = zeros(2^d0t,1);
eta = zeros(2^d0t,1);
psi = zeros(2^d0t,1);
ttimes = zeros(2^d0t, 1);
u = u0;
for t=1:1:2^d0t
    times_0 = tic;

    u_new = mvk(KNm, u, tol, 20, u, 1000); % 1\tI - A/2

    u = mvk(Euler, u, tol, 20, u, 1000);
    
%     norm_rhs = mvk(KNp',u_new,tol,20,u_new,1000);
%     u = u_new;
    for i=1:maxit
        curtime_0 = tic;
        u = dmrg_solve2(KNp, u_new, u, tol, [], [], 1, KNp');
        cur_time = toc(curtime_0);
    
        Mx = mvk(KNp,u,tol,20,u,1000);
        resid_true = norm(Mx-u_new)/norm(u_new);
%         MMx = mvk(KNp',Mx,tol,20,u,1000);
%         resid = norm(MMx - norm_rhs)/norm(norm_rhs);
        
        fprintf('\n\n\t cur_time: %g\n\t true_resid: %3.3e\n\t erank: %g\n', cur_time, resid_true, erank(u));
%        pause(0.5);
        if (resid_true<5*tol) 
            break; 
        end;
    end;

    ttimes(t)=toc(times_0);
    
    angles(t)=acos(dot(u,x_ex)/(norm(u)*norm(x_ex)));
    nrm_u = dot(u,ons);
    tt = zeros(dpx,dpx);
    tt(1,1)=dot(x1.*v1, u);
    tt(2,1)=dot(x2.*v1, u);
    tt(3,1)=dot(x3.*v1, u);
    tt(1,2)=dot(x1.*v2, u);
    tt(2,2)=dot(x2.*v2, u);
    tt(3,2)=dot(x3.*v2, u);
    tt(1,3)=dot(x1.*v3, u);
    tt(2,3)=dot(x2.*v3, u);
    tt(3,3)=dot(x3.*v3, u);
    tt = tt*2/nrm_u; % 2 is dQ^*/dQ
    
    eta(t) = -tt(1,2)/beta;
    psi(t) = -(tt(1,1)-tt(2,2))/(beta^2);
    
    fprintf('\nTime step %d (%3.5e) done. Au/u: %3.3e. Angle(u,x_ex): %3.3e. \n eta: %3.5e. Psi: %3.5e. ttimes: %3.5f\n', t, t*tau, norm(Ax*u)/norm(u), angles(t), eta(t), psi(t), ttimes(t));
    results(t) = norm(Ax*u)/norm(u);
%     pause(1);
%     keyboard;
end;

fprintf('\t d0x=%d\n', d0x);
fprintf('\t d0t=%d\n', d0t);
fprintf('\t eta(T)=%3.15e\n', eta(2^d0t));
fprintf('\t Psi(T)=%3.15e\n', psi(2^d0t));

