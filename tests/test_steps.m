d0t = 7; % quantics dims for t
d0x = 10; % quantics dims for x
dpx = 3; % phys. dims for x

a = -10;
b = 10;
h = (b-a)/(2^d0x+1);

tol = 1e-6;
eps = 1e-9;
maxit = 20;

T = 10;
tau = T/(2^d0t+1);

Ax = tt_matrix(tt_qlaplace_dd(d0x*ones(1,dpx)));
Ax = Ax/(h^2);

% For Fokker-Plank
beta = 1;
ddd = 0.5;
zzz = 0.1;

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
r2 = round(r2, eps);
% eexp = tt_tensor(tt_random(2,r2.d,2));
% for i=1:20
%     eexp = eexp + tt_tensor(tt_random(2,r2.d,2));
%     eexp = funcrs(r2, @(r)exp(-0.5*r/(ddd^2)), 1e-13, eexp, 5);
% end;



v1 = (eexp.*x1)*(zzz/(2*ddd^5)) - x1; % + beta*x2;
v1 = round(v1, eps);
v2 = (eexp.*x2)*(zzz/(2*ddd^5)) - x2; % - x2.*x2.*x2;
v3 = (eexp.*x3)*(zzz/(2*ddd^5)) - x3;
v1 = diag(v1);
v2 = diag(v2);
v3 = diag(v3);

phi = r2*0.5 + 0.5*zzz/(ddd^3)*eexp; % + x2.*x2.*x2.*x2/4;
phi = round(phi, eps);
x_ex = tt_tensor(tt_random(2,r2.d,2));
for i=1:10
    x_ex = x_ex + tt_tensor(tt_random(2,r2.d,2));
    x_ex = funcrs2(phi, @(r)exp(-r), eps, x_ex, 10);
end;

Ax = Ax + Cx1*v1 + Cx2*v2 + Cx3*v3;
Ax = round(Ax, eps);

norm_Ax = norm(Ax*x_ex)/norm(x_ex)

u0 = x_ex;


% Now, generate parts of KN scheme:
KNm = tt_matrix(tt_eye(2,dpx*d0x))/tau - Ax*0.5;
KNm = round(KNm, eps);

KNp = tt_matrix(tt_eye(2,dpx*d0x))/tau + Ax*0.5;
KNp = round(KNp, eps);

u = u0;
for t=1:1:2^d0t
    u_new = mvk(KNm, u, tol, 20, u, 1000); % 1\tI - A/2
    
    norm_rhs = mvk(KNp',u_new,tol,20,u_new,1000);
    u = u_new;
    for i=1:maxit
        tic;
        u = dmrg_solve2(KNp, u_new, u, tol, [], [], 1, KNp');
        cur_time = toc;
    
        Mx = mvk(KNp,u,tol,20,u,1000);
        resid_true = norm(Mx-u_new)/norm(u_new);
        MMx = mvk(KNp',Mx,tol,20,u,1000);
        resid = norm(MMx - norm_rhs)/norm(norm_rhs);
        
        fprintf('\n\n\t cur_time: %g\n\t true_resid: %3.3e\n\t norm.resid: %3.3e\n\t erank: %g\n', cur_time, resid_true, resid, erank(u));
        pause(0.5);
        if (resid_true<5*tol) 
            break; 
        end;
    end;
    
    fprintf('\nTime step %d (%3.5e) done. Dist. to u0: %3.3e\n', t, t*tau, norm(u-u0)/norm(u0));
    pause(1);
end;