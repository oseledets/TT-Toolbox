% d0t = 7; % quantics dims for t
% d0x = 8; % quantics dims for x
dpx = 3; % phys. dims for x

a = -10;
b = 10;
h = (b-a)/(2^d0x+1);

tol = 1e-6;
eps = 1e-8;
maxit = 8;


% start_T = [0, 0.5, 1, 2, 5];
% end_T =   [0.5, 1, 2, 5, 10];

start_T = [0    1.25 2.5  3.75 5.0  6.25 7.5  8.75];
end_T =   [1.25 2.5  3.75 5.0  6.25 7.5  8.75 10.0];


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
eexp = funcrs(r2, @(r)exp(-0.5*r/(ddd^2)), eps, r2, 10);
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
x_ex = funcrs(phi, @(r)exp(-r), eps, phi, 10);

Ax = Ax*0.5 + Cx1*vv1 + Cx2*vv2 + Cx3*vv3;
Ax = round(Ax, eps);

norm_Ax = norm(Ax*x_ex)/norm(x_ex)

% prepare first u0
u0 = x_ex.*x_ex;
u0 = round(u0, eps);

norm_Au0 = norm(Ax*u0)/norm(u0)
% keyboard;


% Now letz init the combined timestep procedure
global_results = cell(1,max(size(start_T)));
norms_Au = cell(1,max(size(start_T)));
for out_t=1:max(size(start_T))
    tau = (end_T(out_t)-start_T(out_t))/(2^d0t);
    
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
    
    Ix3 = kron(Ix, kron(Ix, Ix));
    
    % global matrix
    M = kron(Ix3, Grad_t) + kron(Ax, KN_term);
    M = round(M, eps);
    
    % f1 = sin(pi*(1:1:2^d0x)'/(1+2^d0x));
    % u0 = full_to_qtt(f1, 1e-12);
    % u0 = tt_tensor(u0);
    % u0 = kron(u0,u0);
    
    % u0 = tt_tensor(tt_ones(d0x*dpx, 2));
    
    
    u0_rhs = u0/tau - (Ax*u0)*0.5; % stuff u0 to rhs of KN scheme
    u0_rhs = round(u0_rhs, eps);
    rhs = kron(u0_rhs, e1t);
    
    % norm_rhs = mvk(M',rhs,tol,20,tt_tensor(tt_random(2,rhs.d,2)),1000);
    
    U = tt_random(2, rhs.d, 2);
    
    results = zeros(maxit,6);
    resid_old = 1e15;
    for i=1:maxit
        tic;
        U = dmrg_solve2(M, rhs, U, tol, [], [], 1, [], false);
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
        pause(1);
        
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
    
    global_results{out_t}=results;

    % prepare new start
    last_indices = cell(1, d0t);
    for i=1:d0t
        last_indices{i}=2;
    end;
    u0 = core(U);
    u0 = tt_elem3(u0, last_indices);
    u0 = tt_squeeze(u0);
    u0 = tt_tensor(u0);
    
    norm_Au0 = norm(Ax*u0)/norm(u0)
    norms_Au{out_t} = norm_Au0;
%   keyboard;
end;

    ons = tt_tensor(tt_ones(eexp.d, 2));
    nrm_u = dot(u0,ons);
    tt = zeros(dpx,dpx);
    tt(1,1)=dot(x1.*v1, u0);
    tt(2,1)=dot(x2.*v1, u0);
    tt(3,1)=dot(x3.*v1, u0);
    tt(1,2)=dot(x1.*v2, u0);
    tt(2,2)=dot(x2.*v2, u0);
    tt(3,2)=dot(x3.*v2, u0);
    tt(1,3)=dot(x1.*v3, u0);
    tt(2,3)=dot(x2.*v3, u0);
    tt(3,3)=dot(x3.*v3, u0);
    tt = tt*2/nrm_u; % 2 because we multiplied V by 0.5, buttheads

    eta = -tt(1,2)/beta
    psi = -(tt(1,1)-tt(2,2))/(beta^2)


