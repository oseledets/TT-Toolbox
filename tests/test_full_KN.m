d0t = 12; % quantics dims for t
d0x = 8; % quantics dims for x
dpx = 3; % phys. dims for x
dconf = 4;

a = -10;
b = 10;
h = (b-a)/(2^d0x);

tol = 1e-3;
eps = 1e-8;
maxit = 1;


% start_T = [0, 0.5, 1, 2, 5];
% end_T =   [0.5, 1, 2, 5, 10];

% Trange = [0, 0.5, 2.5, 3.75, 5.0, 6.25, 7.5, 8.75, 10.0];
Trange = [0, 0.5, 1, 2, 5, 10, 20, 50, 80, 100];

% For Fokker-Plank
beta = 0.25;
ddd = 0.5;
zzz = 0.1;

Arouse = spdiags(ones(dconf,1)*[-1,2,-1], [-1,0,1], dconf,dconf);
Arouse = full(Arouse);

% Z = eye(dpx*dconf);

lp1 = tt_matrix(tt_qlaplace_dd(d0x));
lp1 = lp1/(h^2);

Sx = tt_matrix(tt_shf(d0x));
Grad_x = (Sx - Sx')/(2*h);
Grad_x = round(Grad_x, eps);
% Grad_x = tt_matrix(IpaS(d0x,-1));
% Grad_x = Grad_x/h;
Ix = tt_matrix(tt_eye(2, d0x));
% Ix3 = kron(Ix, kron(Ix, Ix));

% Cx1 = kron(Ix, kron(Ix, Grad_x));
% Cx2 = kron(Ix, kron(Grad_x, Ix));
% Cx3 = kron(Grad_x, kron(Ix, Ix));

% Ax = Cx1'*Cx1+Cx2'*Cx2+Cx3'*Cx3;

x = a + (0:1:2^d0x-1)'*h;
eexp1 = exp(-0.5*(x.^2)/(ddd^2));
eexp1 = full_to_qtt(eexp1, eps);
eexp1 = tt_tensor(eexp1);

ttx = a*tt_tensor(tt_ones(d0x,2))+tt_tensor(tt_x(d0x,2))*h;
ex = tt_tensor(tt_ones(d0x,2));
% x1 = kron(kron(ex,ex), ttx);
% x2 = kron(kron(ex,ttx), ex);
% x3 = kron(kron(ttx,ex), ex);

Z = reshape(Z, dconf, dpx, dconf, dpx);


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
                X{i,j}=X{i,j}+Z(i,j,i2,j2)*Xst{i2,j2}; % X = inv(Z')*Xst
            end;
        end;
        X{i,j}=round(X{i,j}, eps);
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
                    curlp=kron(curlp, Ix);
                end;
            end;
        end;
        diaglp{i}=diaglp{i}+curlp;
        diaglp{i}=round(diaglp{i},eps);
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

% Velocities : d\phi / d Qst_i
V = cell(dconf,dpx);
for i=1:dconf
    for j=1:dpx
        V{i,j}=X{i,j};
%         Here we have 0.5*2, two occurences of diag. term
%         V{i,j}=V{i,j} - (zzz/ddd^5)*(diageexp{i}.*Xst{i,j});
%         V{i,j}=round(V{i,j}, eps);
    end;
end;

% Velocities to stuff into equation %%% (Rouse matrix, flow included)
Veq = cell(dconf,dpx);
for i=1:dconf
    for j=1:dpx
        for k=1:dconf            
            Veq{i,j} = Veq{i,j} - 0.25*Arouse(i,k)*V{k,j};            
        end;
        if (j==1)
            Veq{i,j}=Veq{i,j} + beta*X{i,2};
        end;
        Veq{i,j}=round(Veq{i,j}, eps);
    end;
end;

% Stiffness matrix
Ax = [];
for i=1:dconf
    for k=1:dconf
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


% eexp = kron(kron(eexp1, eexp1), eexp1);

% r2 = x1.*x1 + x2.*x2 + x3.*x3;
% r2 = round(r2, eps);
% eexp = funcrs(r2, @(r)exp(-0.5*r/(ddd^2)), eps, r2, 10);
% eexp = funcross(x2.*x2, @(r)exp(-0.5*r/(ddd^2)), eps, x2, 10);



% % v1 = 0.5*(eexp.*x1)*(zzz/(2*ddd^5)) - x1*0.5;
% v1 = 0.5*((eexp.*x1)*(zzz/(2*ddd^5)) - x1);
% v1 = round(v1, eps);
% v2 = 0.5*((eexp.*x2)*(zzz/(2*ddd^5)) - x2); % - x2.*x2.*x2;
% % v3 = 0.5*(eexp.*x3)*(zzz/(2*ddd^5)) - x3*0.5;
% v3 = 0.5*((eexp.*x3)*(zzz/(2*ddd^5)) - x3);

% vv1 = diag(v1 + beta*x2);
% vv2 = diag(v2);
% vv3 = diag(v3);
% 
% phi = r2*0.5 + 0.5*zzz/(ddd^3)*eexp; % + x2.*x2.*x2.*x2/4;
% phi = round(phi, eps);
% x_ex = funcrs(phi, @(r)exp(-r), eps, phi, 10);
% 
% Ax = Ax*0.5 + Cx1*vv1 + Cx2*vv2 + Cx3*vv3;
% Ax = round(Ax, eps);
% 
% norm_Ax = norm(Ax*x_ex)/norm(x_ex)

% prepare first u0
uSN = tt_tensor(full_to_qtt(exp(-0.5*(x.^2)), eps));
% uDD = tt_tensor(full_to_qtt([zeros(2^(d0x-1),1); 1; zeros(2^(d0x-1)-1,1)], eps));
u0=uSN;
for i=2:dpx*dconf
    u0 = kron(u0,uSN);
end;

% u0 = x_ex;
% u0 = round(u0, eps);

norm_Au0 = norm(Ax*u0)/norm(u0)
% keyboard;


% Now letz init the combined timestep procedure
Nt = max(size(Trange))-1;
global_results = cell(Nt+1,1);
norms_Au = zeros(Nt+1,1);
etas = zeros(Nt+1,1);
psis = zeros(Nt+1,1);
for out_t=1:Nt
    tau = (Trange(out_t+1)-Trange(out_t))/(2^d0t);
    
    St = tt_shf(d0t); St = tt_matrix(St); St=St'; % lower shift matrix for gradient
    It = tt_matrix(tt_eye(2,d0t));
    Grad_t = (It-St)/tau; % for gradient
    Grad_t = round(Grad_t, eps);
    
    KN_term = (It+St)*0.5; % Krank-Nikolson term
    KN_term = round(KN_term, eps);
    
    e1t = cell(d0t,1);
    for i=1:d0t
        e1t{i}=[1;0];
    end;
    e1t = tt_tensor(e1t); % first identity vector for t - we need to put u0 into rhs   
    
    % global matrix
    M = kron(tt_matrix(tt_eye(2,d0x*dpx*dconf)), Grad_t) + kron(Ax, KN_term);
    M = round(M, eps);
    
    % f1 = sin(pi*(1:1:2^d0x)'/(1+2^d0x));
    % u0 = full_to_qtt(f1, 1e-12);
    % u0 = tt_tensor(u0);
    % u0 = kron(u0,u0);
    
    % u0 = tt_tensor(tt_ones(d0x*dpx, 2));
        
    u0_rhs = u0/tau - (Ax*u0)*0.5; % stuff u0 into rhs of KN scheme
    u0_rhs = round(u0_rhs, tol);
    rhs = kron(u0_rhs, e1t);
    
    % norm_rhs = mvk(M',rhs,tol,20,tt_tensor(tt_random(2,rhs.d,2)),1000);
    
    U = tt_random(2, rhs.d, 2);
%     U = kron(u0, tt_tensor(tt_ones(d0t,2)));
    
    results = zeros(maxit,6);
    resid_old = 1e15;
    for i=1:maxit
        tic;
        U = dmrg_solve2(M, rhs, tol, 'x0',U, 'nswp', 30, 'verb', 1, 'nrestart', 50, 'max_full_size', 1000);
        cur_time = toc;
        
        if (i==1)
            U_best = U;
        end;
        
        Mx = mvk(M,U,tol,20,tt_tensor(tt_random(2,rhs.d,2)),1000);
        %     Mx = M*x;
        resid_true = norm(Mx-rhs)/norm(rhs);
        %     MMx = mvk(M',Mx,tol,20,tt_tensor(tt_random(2,rhs.d,2)),1000);
        %     resid = norm(MMx - norm_rhs)/norm(norm_rhs);
        resid = resid_true;
        
        if (i>1)
            if (resid_true<results(i-1,2))
                U_best = U;
            else
                U=U_best;
            end;
        end;
        
        fprintf('\n\n\t cur_time: %g\n\t true_resid: %3.3e\n\t norm.resid: %3.3e\n\t erank: %g\n', cur_time, resid_true, resid, erank(U));
        results(i,1)=i;
        results(i,2)=resid_true;
        results(i,3)=resid;
        results(i,4)=erank(U);
        results(i,5)=cur_time;
        pause(0.5);
        
        if (resid_true<5*tol)
            break;
        end;
	if ((resid_true/resid_old>1-1e-2)&&(resid_true<tol*1000))
	    break;
	end;
	resid_old = resid_true;
    end;
    
    results(:,6)=cumsum(results(:,5));
    
    fprintf('sweep\t true resid\t norm. resid \t erank  \t sw. time \t full time\n');
    for i=1:maxit
        fprintf('%d\t %3.3e\t %3.3e\t %3.3f  \t %3.5f \t %3.5f\n', results(i,1),results(i,2),results(i,3),results(i,4),results(i,5),results(i,6));
    end;
    
    global_results{out_t+1}=results;

    % prepare new start
    last_indices = cell(1, d0t);
    for i=1:d0t
        last_indices{i}=2;
    end;
    u0 = core(U);
    u0 = tt_elem3(u0, last_indices);
    u0 = tt_squeeze(u0);
    u0 = tt_tensor(u0);    
    
    ons = tt_tensor(tt_ones(u0.d, 2));
    nrm_u = dot(u0,ons);
    tt = zeros(dpx,dpx);
    for i=1:dconf
        for j=1:dpx
            for k=1:dpx
                tt(j,k)=tt(j,k)-dot(X{i,j}.*V{i,k}, u0);
            end;
        end;
    end;
%     tt(1,1)=dot(x1.*v1, u0);
%     tt(2,1)=dot(x2.*v1, u0);
%     tt(3,1)=dot(x3.*v1, u0);
%     tt(1,2)=dot(x1.*v2, u0);
%     tt(2,2)=dot(x2.*v2, u0);
%     tt(3,2)=dot(x3.*v2, u0);
%     tt(1,3)=dot(x1.*v3, u0);
%     tt(2,3)=dot(x2.*v3, u0);
%     tt(3,3)=dot(x3.*v3, u0);
    tt = tt/nrm_u;

    etas(out_t+1) = -tt(1,2)/beta;
    psis(out_t+1) = -(tt(1,1)-tt(2,2))/(beta^2);
        
    u2 = qtt_to_tt(core(u0), d0x*ones(1,dpx*dconf),0);
    u2 = tt_tensor(u2);
    n=2^d0x/2+1;
    contour(full(u2(:,n,n,n, :,n,n,n, n,n,n,n), 2^d0x*[1,1]));
    
%     Z0 = reshape(Z, dconf*dpx, dconf*dpx);
%     [Z2,x]=normcoords(core(u2), 'x', n*ones(dpx*dconf, 1));
%     Z = Z0*Z2;
    
%     u0f = full(u0, 2^d0x*ones(1,dpx));
    
%     figure(1);
%     mesh(u0f(:,:,1*2^d0x/4));
%     figure(2);
%     mesh(u0f(:,:,2*2^d0x/4));
%     figure(3);
%     mesh(u0f(:,:,3*2^d0x/4));
%     figure(4);
%     mesh(u0f(:,:,4*2^d0x/4));
    
    norms_Au(out_t+1) = norm(Ax*u0)/norm(u0);
    
    fprintf('TimeRange: %d [%g->%g], eta=%3.3e, psi=%3.3e, norm_Au=%3.3e\n', out_t, Trange(out_t), Trange(out_t+1), etas(out_t+1), psis(out_t+1), norms_Au(out_t+1));
    pause(0.5);
%   keyboard;
end;


