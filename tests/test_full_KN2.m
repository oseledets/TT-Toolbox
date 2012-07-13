%The test script for the Crank-Nicolson scheme with global time stepping
times_glob = zeros(9,6);
appr_glob = zeros(9,6);

% for d0t=6:14
% for d0x=6:11
d0t = 15; % quantics dims for t
d0x = 8; % quantics dims for x
dpx = 2; % phys. dims for x

a = 0;
b = 1;
h = (b-a)/(2^d0x+1);


tol = 1e-4;
eps = 1e-8;

% tranges = [0, 0.01, 0.05, 0.1, 0.2, 0.5, 1, 2, 4];
tranges = [0,0.5];


Ax = tt_matrix(tt_qlaplace_dd(d0x*ones(1,dpx)));
Ax = Ax/(h^2);

Ix = tt_matrix(tt_eye(2, dpx*d0x));

x = (a+h:h:b-h)';

z = (0:h:(b-a-2*h))';
z = exp(-z.^2/(4*tranges(end)));
z = tt_reshape(tt_tensor(z), 2*ones(1,d0x), eps);

C1l = tt_matrix(tt_qltrtoepl(core(z)));
C1d = diag(diag(C1l));
C1 = C1l+ C1l' - C1d;
C1 = round(C1, eps)*h;
% C1 = toeplitz(z)*h;
% C1 = tt_matrix(full_to_nd(C1, eps));

C = [];
for i=1:dpx
    C = kron(C,C1);
end;

% u00 = tt_reshape(tt_tensor(sin(pi*x)), 2*ones(1,d0x), eps);
% u00 = tt_ones(2,d0x);
u00 = tt_reshape(tt_tensor(exp(-(x-(a+b)*0.5).^2*16/0.5)), 2*ones(1,d0x), eps);
% u00 = kron(u00,u00);
% u_ex2 = kron(exp(-x.^2*0.5/(0.25+2*tranges(end))), exp(-x.^2*0.5/(0.25+2*tranges(end))))*((4*pi*tranges(end))^(-dpx/2));

x = tt_tensor(x);
x = tt_reshape(x, 2*ones(d0x, 1), eps);
e1 = tt_tensor(tt_ones(2, d0x));
% u0 = zeros(2^d0x,2^d0x);
% u0(2:2^d0x-1, 2:2^d0x-1)=1;
% u00 = tt_reshape(tt_tensor(u00), 2*ones(1,dpx*d0x), eps);
% u_ex2 = tt_reshape(tt_tensor(u_ex2), 2*ones(1,dpx*d0x), eps);
% u0 = tt_tensor(tt_ones(d0x*dpx,2));
% u00 = kron((x-a*e1).*(b*e1-x), e1).*kron(e1, (x-a*e1).*(b*e1-x));

% u00 = ((x-a*e1).*(b*e1-x)/(a*b)).*u00;
% for i=2:gamma
%     u00 = u00.*(x-a*e1).*(b*e1-x)/(a*b);
    u00 = round(u00, eps);
% end;

u00 = kron(u00,u00);

% u_ex = round(C*u00, eps)*((4*pi*tranges(end))^(-dpx/2));
u_ex = u00*exp(-2*pi^2*tranges(end));
% u0 = round(u0, eps);
% keyboard;

spacial_rhs = tt_zeros(2,d0x*dpx);

u0 = u00;
% Now letz init the combined timestep procedure
ttimes = zeros(1, max(size(tranges))-1);
resids = zeros(1, max(size(tranges))-1);
for out_t=1:max(size(tranges))-1
    tau = (tranges(out_t+1)-tranges(out_t))/(2^d0t);
    
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
    eet = tt_tensor(tt_ones(2,d0t));
%     e2t = round(e2t, 1e-12);
    
    
    % global matrix
    M = kron(Ix, Grad_t) + kron(Ax, KN_term);
    M = round(M, eps);
    
    % f1 = sin(pi*(1:1:2^d0x)'/(1+2^d0x));
    % u0 = full_to_qtt(f1, 1e-12);
    % u0 = tt_tensor(u0);
    % u0 = kron(u0,u0);
    
    % u0 = tt_tensor(tt_ones(d0x*dpx, 2));
    
    U = kron(u0, eet);
    
    kickrank = 2;
    trunc_swp = 2;    
    
    tic;
    u0_rhs = u0/tau - (Ax*u0)*0.5; % stuff u0 to rhs of KN scheme
    u0_rhs = round(u0_rhs, eps);
    rhs = kron(u0_rhs, e1t) + kron(spacial_rhs, eet);
    rhs = round(rhs, eps);
   
    % norm_rhs = mvk(M',rhs,tol,20,tt_tensor(tt_random(2,rhs.d,2)),1000);  
%     [U,somedata] = dmrg_solve2(M, rhs, tol, 'x0', U, 'nswp', 20, 'max_full_size', 5000, 'kickrank', 0);
    [U,somedata] = als_adapt_solve(M, rhs, tol, 'x0', U, 'nswp', 20, 'max_full_size', 200, 'kickrank', kickrank, 'trunc_swp', trunc_swp);
    ttimes(out_t) = toc;
        
    somedata_amr_resfactor2{kickrank,trunc_swp} = {ttimes, somedata{4}, somedata{6}};
    
%     Mx = mvk(M,U,tol,20,tt_tensor(tt_random(2,rhs.d,2)),1000);
%     Mx = mvk3(M,U,tol);
    Mx = tt_mvk4(M, U, tol);
    %     Mx = M*x;
    resids(out_t) = norm(Mx-rhs)/norm(rhs);
        

    % prepare new start
    ind = num2cell([1;2]*ones(1,U.d), 1);    
%     last_indices = cell(1, d0t);
    for i=1:d0t
        ind{i+dpx*d0x}=2;
    end;
    u0 = U(ind);
    u0 = tt_reshape(u0, 2*ones(1,dpx*d0x));
%     u0 = core(U);
%     u0 = tt_elem3(u0, last_indices);
%     u0 = tt_squeeze(u0);
%     u0 = tt_tensor(u0);
    
%     appr = norm(u0-u00*exp(-2*pi^2*tranges(out_t+1)))/norm(u00*exp(-2*pi^2*tranges(out_t+1)))
    appr = norm(u0-u_ex)/norm(u_ex)
    
    fprintf('Time range: [%g; %g], solve_time = %g, resid = %3.3e, swps: %d\n', tranges(out_t), tranges(out_t+1), ttimes(out_t), resids(out_t), swps);
%     keyboard;
end;

times_glob(d0t-5, d0x-5)=sum(ttimes);
appr_glob(d0t-5, d0x-5)=appr;

% U2 = tt_reshape(U, [2*ones(1,2*d0x), 2^d0t]);
% ind = cell(2*d0x+1,1); 
% for i=1:2*d0x; ind{i}=':'; end;
% for i=2:2:2^d0t
%     ind{2*d0x+1} = i;
%     mesh(full(U2(ind), 2^d0x*[1,1]))
%     title(sprintf('i=%d', i));
%     drawnow;
%     pause(0.05);
% end;

% end;
% end;