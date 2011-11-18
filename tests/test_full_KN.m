% d0t = 8; % quantics dims for t
d0x = 10; % quantics dims for x
nx = 21;
dpx = 3; % phys. dims for x
dconf = 4;

a = 20; % Domain is [-a,a]^...
b = 10; % For ILangevin scale only!

h = (2*a)/(2^d0x);

tol = 1e-4;
eps = 1e-8;
maxit = 1;


% start_T = [0, 0.5, 1, 2, 5];
% end_T =   [0.5, 1, 2, 5, 10];

% Trange = [0, 0.5, 2.5, 3.75, 5.0, 6.25, 7.5, 8.75, 10.0];
% Trange = [0,0.2, 0.5, 1:1:100];
% d0ts = 10*ones(1,numel(Trange)-1);
Trange = [0, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 80, 100,150,200,300,400];
d0ts =   [ 8,   8,   8,  9, 10, 11, 12, 13, 13, 13,  14, 14, 15, 15];

% For Fokker-Plank
beta = 0.08;
ddd = 0.5;
zzz = 0.1;

Arouse = spdiags(ones(dconf,1)*[-1,2,-1], [-1,0,1], dconf,dconf);
Arouse = full(Arouse);

K = zeros(dpx,dpx);
K(1,2)=beta;
A = kron(eye(dconf), K)-0.25*kron(Arouse, eye(dpx));
A2 = kron(A, eye(dconf*dpx))+kron(eye(dconf*dpx),A);
% f = -2*reshape(eye(dconf*dpx), (dconf*dpx)^2, 1);
f = -reshape(kron(Arouse,eye(dpx)), (dconf*dpx)^2, 1);
P = A2 \ f;
P = reshape(P, dconf*dpx, dconf*dpx);
P = inv(P);

% Z = P^(0.5);
[Z,L]=eig((P+P')/2);
L = diag(L);
L = L.^(-0.5);
L = reshape(L, dpx, dconf);

Z = eye(dpx*dconf);
% [Z,L]=eig(Arouse);
% Z2 = [1,1;-1,1]/sqrt(2);
% Z2 = eye(dpx);
% Z2 = [ -0.607142393387428  -0.794593049398109
%         0.794593049398109  -0.607142393387428];
% Z = kron(Z, Z2);
% 
% lp1 = tt_matrix(tt_qlaplace_dd(d0x));
% lp1 = lp1/(h^2);

% [x,cw]=lgwt(nx, -a, a);
% cw1 = tt_tensor(cw);
% cw = cw1;
% for i=1:dpx*dconf
%     cw = kron(cw, cw1);
% end;
% Grad_x = lagr_diff(x);

% Sx = nd_to_full(tt_shf(d0x));
Sx = tt_matrix(tt_shf(d0x));
Grad_x = (Sx - Sx')/(2*h);
% Grad_x = tt_matrix(Grad_x);
Grad_x = round(Grad_x, eps);
% Grad_x = tt_matrix(IpaS(d0x,-1));
% Grad_x = Grad_x/h;
Ix = tt_matrix(tt_eye(2, d0x));
% Ix = tt_matrix(tt_eye(nx, 1));
% Ix3 = kron(Ix, kron(Ix, Ix));

% Cx1 = kron(Ix, kron(Ix, Grad_x));
% Cx2 = kron(Ix, kron(Grad_x, Ix));
% Cx3 = kron(Grad_x, kron(Ix, Ix));

% Ax = Cx1'*Cx1+Cx2'*Cx2+Cx3'*Cx3;

x = -a + (0:1:2^d0x-1)'*h;
% eexp1 = exp(-0.5*(x.^2)/(ddd^2));
% eexp1 = full_to_qtt(eexp1, eps);
% eexp1 = tt_tensor(eexp1);

% ttx = tt_tensor(x);
ttx = tt_reshape(tt_tensor(x), 2*ones(d0x,1), 1e-10);
% ttx = -a*tt_tensor(tt_ones(d0x,2))+tt_tensor(tt_x(d0x,2))*h;
ex = tt_tensor(tt_ones(d0x,2));
% ex = tt_tensor(tt_ones(1, nx));
% x1 = kron(kron(ex,ex), ttx);
% x2 = kron(kron(ex,ttx), ex);
% x3 = kron(kron(ttx,ex), ex);

% % We need x and y to be in the matrix permutation
% ttx2 = kron2(ttx,ex);
% tty2 = kron2(ex,ttx);
% ex2 = kron2(ex,ex);
% Gradx2 = kron2(Grad_x, Ix);
% Grady2 = kron2(Ix, Grad_x);
% Ix2 = kron2(Ix,Ix);

Gradsst = cell(dconf, dpx);
% Gradient matrices in normal coordinates
for i=1:dconf
    for j=1:dpx
        for i2=1:dconf
            j2=1;
            while (j2<=dpx)
%                 if (i2==i)
%                     if (dpx>1)
%                         if (j==1)&&(j2==1)
%                             Gradsst{i,j}=kron(Gradsst{i,j}, Gradx2);
%                         end;
%                         if (j==2)&&(j2==1)
%                             Gradsst{i,j}=kron(Gradsst{i,j}, Grady2);
%                         end;                        
%                         if (j>2)&&(j2==1)
%                             Gradsst{i,j}=kron(Gradsst{i,j}, Ix2);
%                         end;
%                         if (j>2)&&(j2==j)
%                             Gradsst{i,j}=kron(Gradsst{i,j}, Grad_x);
%                         end;
%                         if (j2>2)&&(j~=j2)
%                             Gradsst{i,j}=kron(Gradsst{i,j}, Ix);
%                         end;
%                     else
%                         if (j2==j)
%                             Gradsst{i,j}=kron(Gradsst{i,j}, Grad_x);
%                         else
%                             Gradsst{i,j}=kron(Gradsst{i,j}, Ix);
%                         end;
%                     end;
%                 else
%                     if (dpx>1)&&(j2==1)
%                         Gradsst{i,j}=kron(Gradsst{i,j}, Ix2);
%                     else
%                         Gradsst{i,j}=kron(Gradsst{i,j}, Ix);
%                     end;
%                 end;
                if (i2==i)&&(j2==j)
                    Gradsst{i,j}=kron(Gradsst{i,j}, Grad_x);
                else
                    Gradsst{i,j}=kron(Gradsst{i,j}, Ix);
                end;
                j2 = j2+1;
            end;
%             if (dpx>1)&&(j2==1)
%                 j2=2;
%             end;
        end;
    end;
end;

% Q_i
Xst = cell(dconf,dpx); % normal coordinates
for i=1:dconf
    for j=1:dpx
        for i2=1:dconf
            j2=1;
            while (j2<=dpx)
%                 if (i2==i)
%                     if (dpx>1)
%                         if (j==1)&&(j2==1)
%                             Xst{i,j}=kron(Xst{i,j}, ttx2);
%                         end;
%                         if (j==2)&&(j2==1)
%                             Xst{i,j}=kron(Xst{i,j}, tty2);
%                         end;                        
%                         if (j>2)&&(j2==1)
%                             Xst{i,j}=kron(Xst{i,j}, ex2);
%                         end;
%                         if (j>2)&&(j2==j)
%                             Xst{i,j}=kron(Xst{i,j}, ttx);
%                         end;
%                         if (j2>2)&&(j~=j2)
%                             Xst{i,j}=kron(Xst{i,j}, ex);
%                         end;
%                     else
%                         if (j2==j)
%                             Xst{i,j}=kron(Xst{i,j}, ttx);
%                         else
%                             Xst{i,j}=kron(Xst{i,j}, ex);
%                         end;
%                     end;
%                 else
%                     if (dpx>1)&&(j2==1)
%                         Xst{i,j}=kron(Xst{i,j}, ex2);
%                     else
%                         Xst{i,j}=kron(Xst{i,j}, ex);
%                     end;
%                 end;
                if (i2==i)&&(j2==j)
                    Xst{i,j}=kron(Xst{i,j}, ttx);
                else
                    Xst{i,j}=kron(Xst{i,j}, ex);
                end;
                j2 = j2+1;
            end;
%             if (dpx>1)&&(j2==1)
%                 j2=2;
%             end;
        end;
    end;
end;

% Laplacen
% diaglp = cell(1,dconf);
% for i=1:dconf
%     for j=1:dpx
%         curlp = [];
%         for j2=1:dpx
%             for i2=1:dconf
%                 if (i2==i)&&(j2==j)
%                     curlp=kron(curlp, lp1);
%                 else
%                     curlp=kron(curlp, Ix);
%                 end;
%             end;
%         end;
%         diaglp{i}=diaglp{i}+curlp;
%         diaglp{i}=round(diaglp{i},eps);
%     end;
% end;

% diagonal repulsions
% diageexp = cell(1,dconf);
% for i=1:dconf
%     for j=1:dpx
%         for i2=1:dconf
%             if (i2==i)
%                 diageexp{i}=kron(diageexp{i}, eexp1);
%             else
%                 diageexp{i}=kron(diageexp{i}, ex);
%             end;
%         end;
%     end;
% end;


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
% uSN = tt_tensor(exp(-0.5*(x.^2)));
% if (dpx>1)
%     uSN2 = funcrs(ttx2.*ttx2+tty2.*tty2, @(x)exp(-0.5*x), eps, ttx2, 20);
%     % uDD = tt_tensor(full_to_qtt([zeros(2^(d0x-1),1); 1; zeros(2^(d0x-1)-1,1)], eps));
%     u0=uSN2;
%     for i=2:dconf
%         u0 = kron(u0, uSN2);
%     end;
%     for i=1:(dpx-2)*dconf
%         u0 = kron(u0,uSN);
%     end;
% else
    u0 = uSN;
    for i=2:dconf*dpx
        u0 = kron(u0,uSN);
    end;
% end;

% u0 = x_ex;
% u0 = round(u0, eps);

% norm_Au0 = norm(Ax*u0)/norm(u0)
% keyboard;

Z0 = zeros(dpx*dconf, dpx*dconf);
% Now letz init the combined timestep procedure
Nt = max(size(Trange))-1;
global_results = cell(Nt+1,1);
norms_Au = zeros(Nt+1,1);
etas = zeros(Nt+1,1);
psis = zeros(Nt+1,1);
Us = cell(Nt,1);
for out_t=1:Nt
    d0t = d0ts(out_t);
    tau = (Trange(out_t+1)-Trange(out_t))/(2^d0t);
    
    
    if (norm(Z(:)-Z0(:))>1e-7) % If Z was updated
        Z = reshape(Z, dpx, dconf, dpx, dconf);
        
        Grads = cell(dconf,dpx); % real gradients
        for j=1:dpx
            for i=1:dconf
                for j2=1:dpx
                    for i2=1:dconf
                        Grads{i,j}=Grads{i,j}+Z(j,i,j2,i2)*Gradsst{i2,j2}; % G=Z*Gst
                    end;
                end;
                Grads{i,j}=round(Grads{i,j}, eps);
            end;
        end;
        
        X = cell(dconf,dpx); % real coordinates
        for j=1:dpx
            for i=1:dconf
                for j2=1:dpx
                    for i2=1:dconf
                        X{i,j}=X{i,j}+Z(j,i,j2,i2)*Xst{i2,j2}; % X = inv(Z')*Xst
                    end;
                end;
                X{i,j}=round(X{i,j}, eps);
            end;
        end;
        
        % Velocities : d\phi / d Qst_i
        V = cell(dconf,dpx);
        for i=1:dconf
            for j=1:dpx
                cx = X{i,j}/b;
                % 7-term Tailor expansion for the Inverse Langevin
                V{i,j}=cx;
%                 V{i,j}=V{i,j}+3/5*(cx.*cx.*cx); % +99/175*(cx.*cx.*cx.*cx.*cx);
%                 V{i,j}=round(V{i,j}, eps);
%                 V{i,j}=V{i,j} + 513/875*(cx.*cx.*cx.*cx.*cx.*cx.*cx);
%                 V{i,j}=round(V{i,j}, eps);
                V{i,j}=V{i,j}*b;
                %             V{i,j}=X{i,j};
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
                %             if (i==k)
                %                 Ax = Ax + 0.25*Arouse(i,k)*diaglp{i};
                %             else
                for j=1:dpx
%                     Ax = Ax + 0.25*Arouse(i,k)*(Grads{i,j}*Grads{k,j}');
                    Ax = Ax + 0.25*Arouse(i,k)*(Grads{i,j}*Grads{k,j}');
                end;
                %             end;
                Ax = round(Ax, eps);
            end;
            for j=1:dpx
                Ax = Ax + Grads{i,j}*diag(Veq{i,j});
            end;
            Ax = round(Ax, eps);
        end;
        
        % Lyapunov function
        Vl = [];
        for i=1:dconf
            for j=1:dpx
                for i2=1:dconf
                    for j2=1:dpx
                        Vl = Vl+(X{i,j}.*X{i2,j2})*P(j+(i-1)*dpx, j2+(i2-1)*dpx);
                        Vl = round(Vl, eps);
                    end;
                end;
            end;
        end;
        u_ex = funcrs2(Vl, @(x)(exp(-x)), tol*0.1, Vl, 25);
%         u_ex = tt_rc(Vl.d, Vl.n, @(ind)(exp(-Vl(ind))), tol, 'x0', u0);
        while (abs(dot(Ax*u_ex, u_ex)/dot(u_ex,u_ex))>1)
            u_ex = funcrs2(Vl, @(x)(exp(-x)), eps, Vl, 25);        
        end;
%         u0 = u_ex;
    end;
        
    Z0 = reshape(Z, dconf*dpx, dconf*dpx);
    
%     St = tt_shf(d0t); St = tt_matrix(St); St=St'; % lower shift matrix for gradient
    Grad_t = IpaS(d0t,-1);
    Grad_t = tt_matrix(Grad_t)/tau;
    It = tt_matrix(tt_eye(2,d0t));
    
%     G2 = kron2(Grad_t, It);
%     iGrad_t = dmrg_solve2(G2, It.tt, 1e-10, 'nswp', 50);
%     iGrad_t = tt_matrix(iGrad_t, 2, 2);
%     Prec = kron(tt_matrix(tt_eye(Ax.n,Ax.d)), iGrad_t);
    
%     It = tt_matrix(tt_eye(2,d0t));
%     Grad_t = (It-St)/tau; % for gradient
    
    KN_term = IpaS(d0t,1);
    KN_term = tt_matrix(KN_term)*0.5; % Krank-Nikolson term
    
    e1t = cell(d0t,1);
    for i=1:d0t
        e1t{i}=[1;0];
    end;
    e1t = tt_tensor(e1t); % first identity vector for t - we need to put u0 into rhs   
    
    % global matrix
    M = kron(tt_matrix(tt_eye(Ax.n,Ax.d)), Grad_t) + kron(Ax, KN_term);
    M = round(M, eps);
    
    % f1 = sin(pi*(1:1:2^d0x)'/(1+2^d0x));
    % u0 = full_to_qtt(f1, 1e-12);
    % u0 = tt_tensor(u0);
    % u0 = kron(u0,u0);
    
    % u0 = tt_tensor(tt_ones(d0x*dpx, 2));
    KNm = tt_matrix(tt_eye(Ax.n,Ax.d))/tau - Ax*0.5;
            
    u0_rhs = mvk3(KNm, u0, tol, 'nswp', 20); % stuff u0 into rhs of KN scheme
%     u0_rhs = u0/tau - (Au0)*0.5; 
%     u0_rhs = round(u0_rhs, tol);
    rhs = kron(u0_rhs, e1t);
    
    % norm_rhs = mvk(M',rhs,tol,20,tt_tensor(tt_random(2,rhs.d,2)),1000);
    
    % Prepare Euler predictor
%     xt = tt_reshape(tt_tensor((1:1:2^d0t)'), 2*ones(d0t,1), eps);
%     et = tt_tensor(tt_ones(d0t,2));
%     xt = tt_matrix(xt, 2,1)*tt_matrix(e1t, 1, 2);
%     et = tt_matrix(et, 2,1)*tt_matrix(e1t, 1, 2);
%     Euler = kron(tt_matrix(tt_eye(Ax.n,Ax.d)), et) - tau*kron(Ax, xt);
%     U = kron(u0, tt_tensor(tt_ones(d0t,2)));
%     U = mvk3(Euler, U, tol);
    
%     U = tt_random(2, rhs.d, 2);
    U = kron(u0, tt_tensor(tt_ones(d0t,2)));
    
    results = zeros(maxit,6);
    resid_old = 1e15;
    for i=1:maxit
        tic;        
        U = dmrg_solve2(M, rhs, tol, 'x0',U, 'nswp', 50, 'verb', 1, 'nrestart', 25, 'max_full_size', 1500, 'min_dpow', 1);
        cur_time = toc;
        
        if (i==1)
            U_best = U;
        end;
        
        Mx = mvk3(M,U,tol*5,'nswp',20); % tt_tensor(tt_random(2,rhs.d,2)),1000
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
    
    Us{out_t}=U;
    
    results(:,6)=cumsum(results(:,5));
    
    fprintf('sweep\t true resid\t norm. resid \t erank  \t sw. time \t full time\n');
    for i=1:maxit
        fprintf('%d\t %3.3e\t %3.3e\t %3.3f  \t %3.5f \t %3.5f\n', results(i,1),results(i,2),results(i,3),results(i,4),results(i,5),results(i,6));
    end;
    
    global_results{out_t+1}=results;

    % prepare new start
    ind = num2cell([1;2]*ones(1,d0x*dpx*dconf), 1);
%     ind = num2cell([1:2^d0x]'*ones(1,dpx*dconf), 1);
%     ind(dpx*dconf+1:dpx*dconf+d0t) = num2cell(2*ones(1,d0t), 1);
    ind(d0x*dpx*dconf+1:d0x*dpx*dconf+d0t) = num2cell(2*ones(1,d0t), 1);
    u0 = U(ind);
%     u0 = tt_reshape(u0, 2^d0x*ones(dpx*dconf, 1));    
    u0 = tt_reshape(u0, 2*ones(d0x*dpx*dconf, 1));    
    
    ons = tt_tensor(tt_ones(u0.d, u0.n));
    nrm_u = dot(u0,ons);
%     nrm_u = dot(u0,cw);
    tt = zeros(dpx,dpx);
    for i=1:dconf
        for j=1:dpx
            for k=1:dpx
                tt(j,k)=tt(j,k)-dot(X{i,j}.*V{i,k}, u0);
%                 tt(j,k)=tt(j,k)-dot(cw.*X{i,j}.*V{i,k}, u0);
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

%     etas(out_t+1)=tt(1,1);
    etas(out_t+1) = -tt(1,2)/beta;
    psis(out_t+1) = -(tt(1,1)-tt(2,2))/(beta^2);
        
    u2 = tt_reshape(u0, 2^d0x*ones(dpx*dconf, 1));
%     u2 = u0;
    n=2^d0x/2+1;
    plotind = 1:max(2^(d0x-7),1):2^d0x;    
    ind = num2cell(n*ones(1, dpx*dconf), 1);
    ind(1)={plotind};
    ind(2)={plotind};
    figure(1);
    contour(full(u2(ind), numel(plotind)*[1,1]));
%     ind(3)={plotind};
%     ind(4)={plotind};
%     ind(1)={n};
%     ind(2)={n};
%     figure(2);
%     contour(full(u2(ind), numel(plotind)*[1,1]));    
%     figure(1);
%     contour(full(u2(plotind,plotind,n,n,n,n,n), 2^7*[1,1]));   
%     figure(2);
%     contour(full(u2(plotind,n,n,plotind), 2^7*[1,1]));   
    
    
%     [Z2,x]=normcoords(core(u2), 'x', n*ones(dpx*dconf, 1));
%     Z = Z0*Z2;
    
%     gridold = cell(dpx*dconf,1);
%     for i=1:dpx*dconf
%         gridold{i}=qtt_to_full(core(ttx));
%     end;
%     Z2 = reshape(Z2, dconf, dpx, dconf, dpx);
%     X = cell(dconf,dpx); % new coordinates
%     for j=1:dpx
%         for i=1:dconf
%             for j2=1:dpx
%                 for i2=1:dconf
%                     X{i,j}=X{i,j}+Z2(i,j,i2,j2)*Xst{i2,j2}; % X = inv(Z')*Xst
%                 end;
%             end;
%             X{i,j}=round(X{i,j}, eps);
%             X{i,j}=core(X{i,j});
%         end;
%     end;    
%     u2 = core(u2);
%     u0 = tt_rc(d0x*dpx*dconf, 2, @(ind)interp_u(ind,u2,X,gridold,dpx,dconf), tol);
    
%     u0f = full(u0, 2^d0x*ones(1,dpx));
    
%     figure(1);
%     mesh(u0f(:,:,1*2^d0x/4));
%     figure(2);
%     mesh(u0f(:,:,2*2^d0x/4));
%     figure(3);
%     mesh(u0f(:,:,3*2^d0x/4));
%     figure(4);
%     mesh(u0f(:,:,4*2^d0x/4));
    
    Au0 = mvk3(tt_matrix(tt_eye(Ax.n,Ax.d))+Ax, u0, tol, 'nswp', 20);
    norms_Au(out_t+1) = norm(Au0-u0)/norm(u0);
    
    fprintf('TimeRange: %d [%g->%g], eta=%3.3e, psi=%3.3e, norm_Au=%3.3e\n', out_t, Trange(out_t), Trange(out_t+1), etas(out_t+1), psis(out_t+1), norms_Au(out_t+1));
    pause(0.5);
%   keyboard;
end;




% 
% for t=1:Nt
% d0t = d0ts(out_t);
% ind = num2cell([1;2]*ones(1,d0x*dpx*dconf), 1);
% ind(d0x*dpx*dconf+1:d0x*dpx*dconf+d0t) = num2cell(2*ones(1,d0t), 1);
% u0 = Us{t}(ind);
% u0 = tt_reshape(u0, 2*ones(d0x*dpx*dconf, 1));
% fprintf('%3.5e \t %g\n', dot(Ax*u0,u0)/dot(u0,u0), erank(round(u0, tol)));
% end;