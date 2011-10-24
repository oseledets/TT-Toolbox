d0t = 9; % quantics dims for t
d0x = 9; % quantics dims for x
dpx = 3; % phys. dims for x
dconf = 4;

a = -20;
b = 20;
% a=0; b=1;
h = (b-a)/(2^d0x+1);

tol = 1e-5;
eps = 1e-8;
maxit = 5;

T = 20;
tau = T/(2^d0t);

Arouse = spdiags(ones(dconf,1)*[-1,2,-1], [-1,0,1], dconf,dconf);
Arouse = full(Arouse);

% Z = eye(dconf*dpx);
% for i=1:dconf
%     Z(:,i)=zeros(dconf*dpx,1);
% %     for j=1:dpx
%         Z((1-1)*dconf+i,i)=1;
%         Z((2-1)*dconf+i,i)=1;
% %     end;
% end;
% [Z,rv]=qr(Z);
% Z = eye(dconf*dpx);
% Z = cell(dpx,1);
% iZ = cell(dpx,1);
% Z{1}=eye(dconf);
% Z{1}(:,1)=[1.974460611288290; 2.193029490616621; 2.000000000000000];
% Z{1}=[-0.491545316295007,  -0.703122019335240,  -0.513811860464874;
%    0.709191616987219,   0.019215831067080,  -0.704753859323562;
%   -0.505401278612842,   0.710809522841623,  -0.489202790071938];
% Z{1}=[ones(dconf,1), [eye(dconf-1); zeros(1,dconf-1)]];
% [Z{1},rv]=qr(Z{1});
% G2 = [1,1;-1,1]/sqrt(2);
% G3 = [0.406473364091040, 0.913662631546520; -0.913662631546520, 0.406473364091040];
% Z{1}=[eye(1),zeros(1,3);zeros(2,1),G3,zeros(2,1);zeros(1,3),eye(1)]*[G2,zeros(2,2);zeros(2,2),eye(2)]; % *[eye(2),zeros(2,2);zeros(2,2),G2];
% iZ{1} = inv(Z{1});
% for i=2:dpx
%     Z{i} = eye(dconf);
%     iZ{i} = eye(dconf);
% end;
% [V,L]=eig(Arouse);
% Z = V*L^(-0.5)*2;
% iZ = inv(Z);

iZ = inv(Z');
Z = reshape(Z, dconf, dpx, dconf, dpx);
iZ = reshape(iZ, dconf, dpx, dconf, dpx);

lp1 = tt_matrix(tt_qlaplace_dd(d0x));
lp1 = lp1/(h^2);

% For Fokker-Plank
beta = 0.5;
ddd = 1;
zzz = 0.1;


Sx = tt_matrix(tt_shf(d0x));
Grad_x = (Sx - Sx')/(2*h);
Grad_x = round(Grad_x, eps);
% Grad_x = tt_matrix(IpaS(d0x,-1));
% Grad_x = Grad_x/h;
Ix = tt_matrix(tt_eye(2, d0x));

Cx1 = kron(Ix, kron(Ix, Grad_x));
Cx2 = kron(Ix, kron(Grad_x, Ix));
Cx3 = kron(Grad_x, kron(Ix, Ix));

x = a + (0:1:2^d0x-1)'*h;
eexp1 = exp(-0.5*(x.^2)/(ddd^2));
eexp1 = full_to_qtt(eexp1, eps);
eexp1 = tt_tensor(eexp1);

ttx = tt_tensor(full_to_qtt(x, eps));
% ttx = a*tt_tensor(tt_ones(d0x,2))+tt_tensor(tt_x(d0x,2))*h;
ex = tt_tensor(tt_ones(d0x,2));
x1 = kron(kron(ex,ex), ttx);
x2 = kron(kron(ex,ttx), ex);
x3 = kron(kron(ttx,ex), ex);

ex3 = tt_tensor(tt_ones(d0x*dpx,2));
% eexp = kron(kron(eexp1, eexp1), eexp1);

% r2 = x1.*x1 + x2.*x2 + x3.*x3;
% r2 = round(r2, eps);
% eexp = funcrs(r2, @(r)exp(-0.5*r/(ddd^2)), eps, r2, 10);
% eexp = funcross(x2.*x2, @(r)exp(-0.5*r/(ddd^2)), eps, x2, 10);

Gradsst = cell(dconf, dpx);
% Gradient matrices in normal coordinates
for j=1:dpx
    for i=1:dconf
        for j2=1:dpx
            for i2=1:dconf
                if (i2==i)&&(j2==j)
                    Gradsst{i,j}=kron(Gradsst{i,j}, Grad_x);
                else
                    Gradsst{i,j}=kron(Gradsst{i,j}, Ix);
                end;
            end;
        end;
    end;
end;
Grads = cell(dconf,dpx); % real gradients
for j=1:dpx
    for i=1:dconf
        for j2=1:dpx
            for i2=1:dconf
                Grads{i,j}=Grads{i,j}+Z(i,j,i2,j2)*Gradsst{i2,j2}; % G=Z*Gst
            end;
        end;
        Grads{i,j}=round(Grads{i,j}, eps);
    end;
end;

% Q_i
Xst = cell(dconf,dpx); % normal coordinates
for j=1:dpx
    for i=1:dconf
        for j2=1:dpx
            for i2=1:dconf
                if (i2==i)&&(j2==j)
                    Xst{i,j}=kron(Xst{i,j}, ttx);
                else
                    Xst{i,j}=kron(Xst{i,j}, ex);
                end;
            end;
        end;
    end;
end;
X = cell(dconf,dpx); % real coordinates
for j=1:dpx
    for i=1:dconf
        for j2=1:dpx
            for i2=1:dconf
                X{i,j}=X{i,j}+iZ(i,j,i2,j2)*Xst{i2,j2}; % X = inv(Z')*Xst
            end;
        end;
        X{i,j}=round(X{i,j}, eps);
    end;
end;

% diagonal repulsions
diageexp = cell(1,dconf);
for i=1:dconf
    for j=1:dpx
        for i2=1:dconf
            if (i2==i)
                diageexp{i}=kron(diageexp{i}, eexp1);
            else
                diageexp{i}=kron(diageexp{i}, ex);
            end;
        end;
    end;
end;

% Laplacen
diaglp = cell(1,dconf);
for i=1:dconf
    for j=1:dpx
        curlp = [];
        for j2=1:dpx
            for i2=1:dconf
                if (i2==i)&&(j2==j)
                    curlp=kron(curlp, lp1);
                else
                    curlp=kron(curlp, tt_matrix(tt_eye(2, d0x)));
                end;
            end;
        end;
        diaglp{i}=diaglp{i}+curlp;
        diaglp{i}=round(diaglp{i},eps);
    end;
end;

% Velocities : d\phi / d Qst_i
V = cell(dconf,dpx);
% Vfunc = @(x)(x - (zzz/ddd^5)*exp(-0.5*x.^2/(ddd^2)).*x);
for i=1:dconf
%     for k=1:dconf
%         r2 = X{1,k}.*X{1,k}+X{2,k}.*X{2,k}+X{3,k}.*X{3,k};
%         r2 = round(r2, eps);
%         Vfunc = @(x)(1+0-(zzz/ddd^5)*exp(-0.5*x/(ddd^2)));
%         curV = funcrs(r2, Vfunc, tol, r2, 10);
%         for j=1:dpx
% %             Vfunc = @(ind)(X{j,k}(ind)-(zzz/ddd^5)*exp(-0.5*r2(ind)/(ddd^2)).*X{j,k}(ind));
% %             curV = tt_rc(d0x*dpx*dconf, 2*ones(d0x*dpx*dconf,1), Vfunc, eps, 'nswp', 10);
% %             curV = funcrs(X{j,k}, Vfunc, eps, X{j,k}, 10);
%             V{j,i}=V{j,i}+iZ(k,i)*(curV-0*tt_tensor(tt_ones(d0x*dpx*dconf,2))).*X{j,k};
%             V{j,i}=round(V{j,i}, eps);
%         end;  

        for j=1:dpx
            V{i,j}=X{i,j}*2;
%             Here we have 0.5*2, two occurences of diag. term
%             V{i,j}=V{i,j} - (zzz/ddd^5)*(diageexp{i}.*Xst{i,j});
%             V{i,j}=round(V{i,j}, eps);
        end;
%     end;
end;
% Velocities to stuff into equation %%% (Rouse matrix, flow included)
Veq = cell(dconf,dpx);
for i=1:dconf
    for j=1:dpx
%         Veq{j,i}=-V{j,i};
        Veq{i,j}=0*tt_tensor(tt_ones(dpx*dconf*d0x,2));
        for k=1:dconf            
            Veq{i,j} = Veq{i,j} - 0.25*Arouse(i,k)*V{k,j};            
            Veq{i,j}=tt_tensor(tt_compr2(core(Veq{i,j}), eps));
        end;
        if (j==1)
            Veq{i,j}=Veq{i,j} + beta*X{i,2};
            Veq{i,j}=round(Veq{i,j}, eps);
        end;
    end;
end;

% % v1 = 0.5*(eexp.*x1)*(zzz/(2*ddd^5)) - x1*0.5;
% v1 = 0.5*((eexp.*x1)*(zzz/(2*ddd^5)) - x1);
% v1 = round(v1, eps);
% v2 = 0.5*((eexp.*x2)*(zzz/(2*ddd^5)) - x2); % - x2.*x2.*x2;
% % v3 = 0.5*(eexp.*x3)*(zzz/(2*ddd^5)) - x3*0.5;
% v3 = 0.5*((eexp.*x3)*(zzz/(2*ddd^5)) - x3);

% vv1 = diag(v1 + beta*x2);
% vv2 = diag(v2);
% vv3 = diag(v3);

% Phi
% phi = 0*tt_tensor(tt_ones(dpx*dconf*d0x,2));
% for i=1:dconf
%     r2 = (X{1,i}.*X{1,i}+X{2,i}.*X{2,i}+X{3,i}.*X{3,i});
%     r2 = round(r2, eps);
%     phi = phi+funcrs(r2, @(x)(0.5*x+zzz/(ddd^3)*exp(-0.5*x/(ddd^2))), eps, r2, 10);
% %     phi = phi + 0.5*(X{1,i}.*X{1,i}+X{2,i}.*X{2,i}+X{3,i}.*X{3,i});
%     % Here we have 0.5*2, two occurences of diag. term
% %     phi = phi + zzz/(ddd^3)*diageexp{i};
%     phi = round(phi, eps);
% end;

% phi = r2*0.5 + 0.5*zzz/(ddd^3)*eexp; % + x2.*x2.*x2.*x2/4;
% phi = round(phi, eps);
% x_ex = funcrs(phi, @(r)exp(-r), eps, phi, 20);
% x_ex = tt_rc(d0x*dpx*dconf,2,@(ind)exp(-phi(ind)),eps);

% Stiffness matrix
Ax = 0*tt_matrix(tt_eye(2,d0x*dpx*dconf));
for i=1:dconf
    for k=1:dconf
%         Ax = Ax + diaglp{i};
        if (i==k)
            Ax = Ax + 0.25*Arouse(i,k)*diaglp{i};
        else
            Ax = Ax + 0.25*Arouse(i,k)*(Grads{i,1}*Grads{k,1}'+Grads{i,2}*Grads{k,2}'+Grads{i,3}*Grads{k,3}');
        end;
        Ax = round(Ax, eps);
    end;    
    Ax = Ax + Grads{i,1}*diag(Veq{i,1}) + Grads{i,2}*diag(Veq{i,2}) + Grads{i,3}*diag(Veq{i,3});
    Ax = round(Ax, eps);
end;

% Ax2 = lp*0.5 + Cx1*diag(Veq{3}) + Cx2*diag(Veq{2}) + Cx3*diag(Veq{1});
% Ax2 = round(Ax2, eps);

% norm_Ax = norm(Ax*x_ex)/norm(x_ex)
% norm_Ax2 = norm(Ax2*x_ex)/norm(x_ex)

uSN = tt_tensor(full_to_qtt(exp(-0.5*(x.^2)), eps));
% uDD = tt_tensor(full_to_qtt([zeros(2^(d0x-1),1); 1; zeros(2^(d0x-1)-1,1)], eps));
u0=uSN;
for i=2:dpx*dconf
    u0 = kron(u0,uSN);
end;
% u0 = round(u0, tol);
% u0 = x_ex;

% f1 = sin((a+(1:1:2^d0x)*h)'*pi);
% f1 = 

% Now, generate parts of KN scheme:
KNm = tt_matrix(tt_eye(2,dpx*d0x*dconf))/tau - Ax*0.5;
KNm = round(KNm, eps);
KNm = core(KNm);

KNp = tt_matrix(tt_eye(2,dpx*d0x*dconf))/tau + Ax*0.5;
KNp = round(KNp, eps);
KNp = core(KNp);

Euler = tt_matrix(tt_eye(2,dpx*d0x*dconf)) - Ax*tau;
Euler = round(Euler, eps);
Euler = core(Euler);

ons = tt_tensor(tt_ones(dpx*d0x*dconf, 2));

results = zeros(2^d0t,1);
% angles = zeros(2^d0t,1);
eta = zeros(2^d0t,1);
psi = zeros(2^d0t,1);
ttimes = zeros(2^d0t, 1);
u = u0;
for t=1:1:2^d0t
    times_0 = tic;

%     u_new = mvk(KNm, u, tol, 20, u, 1000); % 1\tI - A/2
    u_new = tt_mvk3(KNm, core(u), tol, 'verb', 1, 'nswp', 25, 'y0', core(u), 'kickrank', 1);
    u_new = tt_tensor(u_new);

%     u = mvk(Euler, u, tol, 20, u, 1000);
%     u = tt_mvk3(Euler, core(u), tol, 'verb', 1, 'nswp', 25, 'y0', core(u));
    
%     norm_rhs = mvk(KNp',u_new,tol,20,u_new,1000);
%     u = u_new;
    for i=1:maxit
        curtime_0 = tic;
        u = dmrg_solve2(KNp, u_new, tol, 'x0', u, 'nswp', 10, 'verb', 1, 'max_full_size', 1000, 'min_dpow', 1, 'kickrank', 1);
%         u = tt_gmres(KNp, core(u_new), tol, 10,20, tol, tol, [], [], [], core(u));
%         u = tt_tensor(u);
%         u = tt_tensor(u);
        cur_time = toc(curtime_0);
    
%         Mx = mvk(KNp,u,tol,20,u,1000);
        Mx = tt_mvk3(KNp, core(u), tol, 'verb', 1, 'nswp', 10, 'y0', core(u));
        Mx = tt_tensor(Mx);
        resid_true = norm(Mx-u_new)/norm(u_new);
%         resid_true = NaN;
%         MMx = mvk(KNp',Mx,tol,20,u,1000);
%         resid = norm(MMx - norm_rhs)/norm(norm_rhs);
        
        fprintf('\n\n\t cur_time: %g\n\t true_resid: %3.3e\n\t erank: %g\n', cur_time, resid_true, erank(u));
%        pause(0.5);
        if (resid_true<5*tol) 
            break; 
        end;
    end;

    ttimes(t)=toc(times_0);
    
%     angles(t)=acos(dot(u,x_ex)/(norm(u)*norm(x_ex)));
    nrm_u = dot(u,ons);
    % Kramer funcitonal
    tt = zeros(dpx,dpx);
    for i=1:dconf 
        for j=1:dpx
            for k=1:dpx
%                 curV=[];
%                 for m=1:dconf
%                     curV = curV+Z(m,i)*V{k,m};
%                 end;
                tt(j,k)=tt(j,k)-dot(X{i,j}.*V{i,k}, u);
            end;
        end;
%         tt(1,1)=dot(x1.*v1, u);
%         tt(2,1)=dot(x2.*v1, u);
%         tt(3,1)=dot(x3.*v1, u);
%         tt(1,2)=dot(x1.*v2, u);
%         tt(2,2)=dot(x2.*v2, u);
%         tt(3,2)=dot(x3.*v2, u);
%         tt(1,3)=dot(x1.*v3, u);
%         tt(2,3)=dot(x2.*v3, u);
%         tt(3,3)=dot(x3.*v3, u);
    end;
    tt = tt/nrm_u; % % 2 because we multiplied V by 0.5, buttheads
    
    eta(t) = -tt(1,2)/beta;
    psi(t) = -(tt(1,1)-tt(2,2))/(beta^2);
    
%     figure(1);
%     plot((1:t)*tau, eta(1:t));
%     figure(2);
%     plot((1:t)*tau, psi(1:t));    
    
%     u0f = core(u);
%     u0f = qtt_to_tt(u0f, d0x*ones(1,dpx), 0);
%     u0f{3} = u0f{3}(2*2^d0x/4, :, :);
%     u0f = tt_to_full(u0f);
%     figure(1);
%     mesh(u0f(:,:,1*2^d0x/4));
%     figure(3);
%     mesh(u0f);
%     figure(3);
%     mesh(u0f(:,:,3*2^d0x/4));
%     figure(4);
%     mesh(u0f(:,:,4*2^d0x/4));    

%     results(t) = norm(Ax*u)/norm(u);

    u2 = qtt_to_tt(core(u), d0x*ones(1,dpx*dconf),0);
    u2 = tt_tensor(u2);
%     ind = cell(dpx*dconf-2, 1);
%     for i=1:dpx*dconf-2
%         ind{i}=2^(d0x-1);
%     end;
%     u2 = tt_elem3(u2, ind);
    figure(1)
    contour(full(u2(2^(d0x-1),2^(d0x-1),2^(d0x-1),2^(d0x-1), 2^(d0x-1),2^(d0x-1),2^(d0x-1),2^(d0x-1), 2^(d0x-1),2^(d0x-1),:,:), 2^d0x*[1,1]));
%     view(2);
%     figure(2)
%     contour(full(u2(:,2^(d0x-1),2^(d0x-1),2^(d0x-1), :,2^(d0x-1),2^(d0x-1),2^(d0x-1), 2^(d0x-1),2^(d0x-1),2^(d0x-1),2^(d0x-1)), 2^d0x*[1,1]));
%     view(2);
    
    
    results(t)=NaN;
%     fprintf('\nTime step %d (%3.5e) done. Au/u: %3.3e. Angle(u,x_ex): %3.3e. \n eta: %3.5e. Psi: %3.5e. ttimes: %3.5f\n', t, t*tau, norm(Ax*u)/norm(u), angles(t), eta(t), psi(t), ttimes(t));
    fprintf('\nTime step %d (%3.5e) done. Au/u: %3.3e. \n eta: %3.5e. Psi: %3.5e. ttimes: %3.5f\n', t, t*tau, results(t), eta(t), psi(t), ttimes(t));   
    pause(0.5);
%     keyboard;
end;

fprintf('\t d0x=%d\n', d0x);
fprintf('\t d0t=%d\n', d0t);
fprintf('\t eta(T)=%3.15e\n', eta(2^d0t));
fprintf('\t Psi(T)=%3.15e\n', psi(2^d0t));

