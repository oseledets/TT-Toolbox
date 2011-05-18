d0t = 7;
d0x = 6;
nn = 2;

a = -10;
b = 10;
% a=0; b=1;
h = (b-a)/(2^d0x+1);

tol = 1e-6;
eps = 1e-8;
maxit = 20;

T = 50;
tau = T/(2^d0t+1);

% For Fokker-Plank
beta = 1;
ddd = 0.5;
zzz = 0.1;

ttx = a*tt_tensor(tt_ones(d0x,2))+tt_tensor(tt_x(d0x,2))*h;
ex = tt_tensor(tt_ones(d0x,2));
ex3 = tt_tensor(tt_ones(d0x*3,2));
Q = cell(nn-1,3);
for i=1:nn-1
    Q{i,1}=kron(kron(ttx,ex), ex);
    Q{i,2}=kron(kron(ex,ttx), ex);
    Q{i,3}=kron(kron(ex,ex), ttx);
    for j=1:i-1
        Q{i,1}=kron(ex3, Q{i,1});
        Q{i,2}=kron(ex3, Q{i,2});
        Q{i,3}=kron(ex3, Q{i,3});
    end;
    for j=i+1:nn-1
        Q{i,1}=kron(Q{i,1}, ex3);
        Q{i,2}=kron(Q{i,2}, ex3);
        Q{i,3}=kron(Q{i,3}, ex3);
    end;    
end;
return
eexp = cell(nn,nn);
for mu=1:nn
    for nu=1:nn
        if (mu~=nu)
            curx = 0*tt_tensor(tt_ones(Q{1,1}.d, 2));
            cury = 0*tt_tensor(tt_ones(Q{1,1}.d, 2));
            curz = 0*tt_tensor(tt_ones(Q{1,1}.d, 2));
            for i=min(mu,nu):(max(mu,nu)-1)
                curx = curx + Q{i,1};
                cury = cury + Q{i,2};
                curz = curz + Q{i,3};
            end;
            r2 = curx.*curx + cury.*cury + curz.*curz;
            r2 = round(r2, eps);
            eexp{mu,nu} = funcross(r2, @(r)exp(-0.5*r/(ddd^2)), eps, r2, 10);
        end;        
    end;
end;

V = cell(nn-1,3);

for i=1:nn-1
    V{i,1}=-Q{i,1};
    V{i,2}=-Q{i,2};
    V{i,3}=-Q{i,3};
    
    for mu=i:nn
        if (mu~=i)
            curx = 0*tt_tensor(tt_ones(Q{1,1}.d, 2));
            cury = 0*tt_tensor(tt_ones(Q{1,1}.d, 2));
            curz = 0*tt_tensor(tt_ones(Q{1,1}.d, 2));
            for j=i:mu-1
                curx = curx + Q{j,1};
                cury = cury + Q{j,2};
                curz = curz + Q{j,3};
            end;
            V{i,1}=V{i,1}+((eexp{mu,i}+eexp{i,mu}).*curx)*(zzz/(2*ddd^5));
            V{i,2}=V{i,2}+((eexp{mu,i}+eexp{i,mu}).*cury)*(zzz/(2*ddd^5));
            V{i,3}=V{i,3}+((eexp{mu,i}+eexp{i,mu}).*curz)*(zzz/(2*ddd^5));
        end;
    end;
    
    V{i,1}=V{i,1}/(4*(sin(pi*i/(2*nn)))^2);
    V{i,1}=round(V{i,1}, eps);
    V{i,2}=V{i,2}/(4*(sin(pi*i/(2*nn)))^2);
    V{i,2}=round(V{i,2}, eps);
    V{i,3}=V{i,3}/(4*(sin(pi*i/(2*nn)))^2);
    V{i,3}=round(V{i,3}, eps);
end;

Sx = tt_matrix(tt_shf(d0x));
Grad_x = (Sx - Sx')/(2*h);
Grad_x = round(Grad_x, eps);
Ix = tt_matrix(tt_eye(2, d0x));

Cx1 = kron(Grad_x, kron(Ix, Ix));
Cx2 = kron(Ix, kron(Grad_x, Ix));
Cx3 = kron(Ix, kron(Ix, Grad_x));

Ix3 = tt_matrix(tt_eye(2, d0x*3));
Cx = cell(nn-1,3);
for i=1:nn-1
    Cx{i,1}=Cx1;
    Cx{i,2}=Cx2;
    Cx{i,3}=Cx3;
    for j=1:i-1
        Cx{i,1}=kron(Ix3, Cx{i,1});
        Cx{i,2}=kron(Ix3, Cx{i,2});
        Cx{i,3}=kron(Ix3, Cx{i,3});
    end;
    for j=i+1:nn-1
        Cx{i,1}=kron(Cx{i,1}, Ix3);
        Cx{i,2}=kron(Cx{i,2}, Ix3);
        Cx{i,3}=kron(Cx{i,3}, Ix3);
    end;        
end;

lp = 0*tt_matrix(tt_eye(2, d0x*3*(nn-1)));
for i=1:nn-1
    curlp = tt_matrix(tt_qlaplace_dd(d0x*ones(1,3)));
    curlp=curlp/(h^2);
    curlp = curlp/(4*(sin(pi*i/(2*nn)))^2);
    for j=1:i-1
        curlp=kron(Ix3, curlp);
    end;
    for j=i+1:nn-1
        curlp=kron(curlp, Ix3);
    end;            
    lp = lp+curlp;
end;

Ax = lp;
for i=1:nn-1
    for j=1:3
        Ax = Ax + Cx{i,j}*diag(V{i,j});
    end;
    Ax = Ax+Cx{i,1}*diag(beta*Q{i,2});
    Ax = round(Ax, eps);    
end;

u0 = eexp{1,2};

% Now, generate parts of KN scheme:
KNm = tt_matrix(tt_eye(2,3*d0x*(nn-1)))/tau - Ax*0.5;
KNm = round(KNm, eps);

KNp = tt_matrix(tt_eye(2,3*d0x*(nn-1)))/tau + Ax*0.5;
KNp = round(KNp, eps);

Euler = tt_matrix(tt_eye(2,3*d0x*(nn-1))) - Ax*tau;
Euler = round(Euler, eps);

ons = tt_tensor(tt_ones(u0.d, 2));

Au = zeros(2^d0t,1);
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
    
    nrm_u = dot(u,ons);
    tt = zeros(3,3);
    for i=1:nn-1
        tt(1,1)=tt(1,1)+dot(Q{i,1}.*V{i,1}, u);
        tt(2,1)=tt(2,1)+dot(Q{i,2}.*V{i,1}, u);
        tt(3,1)=tt(3,1)+dot(Q{i,3}.*V{i,1}, u);
        tt(1,2)=tt(1,2)+dot(Q{i,1}.*V{i,2}, u);
        tt(2,2)=tt(2,2)+dot(Q{i,2}.*V{i,2}, u);
        tt(3,2)=tt(3,2)+dot(Q{i,3}.*V{i,2}, u);
        tt(1,3)=tt(1,3)+dot(Q{i,1}.*V{i,3}, u);
        tt(2,3)=tt(2,3)+dot(Q{i,2}.*V{i,3}, u);
        tt(3,3)=tt(3,3)+dot(Q{i,3}.*V{i,3}, u);
    end;
    tt = tt*2/nrm_u; % 2 is dQ^*/dQ
    
    eta(t) = -tt(1,2)/beta;
    psi(t) = -(tt(1,1)-tt(2,2))/(beta^2);
    
    Au(t) = norm(Ax*u)/norm(u);    
    fprintf('\n\tTime step %d (%3.5e) done. Au/u: %3.3e. ttimes: %3.5f \n\t eta: %3.5e. Psi: %3.5e\n', t, t*tau, Au(t), ttimes(t), eta(t), psi(t));
    pause(1);
%     keyboard;
end;
