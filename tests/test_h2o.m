% Chemical master equation for 2H + O <-> H2O

L = 7;  % QTT dimension in state variables

Trange = 0:0.5:5;
d0ts = 12*ones(1,numel(Trange)-1);

tol = 1e-7;
eps = 1e-12;

% initial distribution
l = [13,15,0]; % H, O, H2O

% One shift z=-1 (for O)
S1 = tt_shf(L);
S2 = round(S1*S1, eps); %  double shift z=-2 (for H)
% z=1 for H2O will be S1'

I = tt_eye(2,3*L);
e1 = tt_tensor(num2cell([1;0]*ones(1,L), 1));
e2 = zeros(2^L,1); e2(2)=1; e2 = tt_tensor(reshape(e2, 2*ones(1,L)),eps);

x = tt_x(2,L);
o = tt_ones(2,L);

w1 = 2*mtkron(o-e1-e2, o-e1, o);
w2 = 1*mtkron(o, o, o-e1);
% w = mtkron(x, x, o);

% CME matrix
A = -(mtkron(S2,S1,S2')*diag(w1) - diag(w1)); % forward
A = A - (mtkron(S2',S1',S2)*diag(w2) - diag(w2)); % backward
A = round(A, eps);

% initial state
e1 = zeros(2^L,1); e1(l(1)+1)=1;
u0 = tt_tensor(e1);
e1 = zeros(2^L,1); e1(l(2)+1)=1;
u0 = tkron(u0, tt_tensor(e1));
e1 = zeros(2^L,1); e1(l(3)+1)=1;
u0 = tkron(u0, tt_tensor(e1));
% make QTT
u0 = tt_reshape(u0, 2*ones(1,3*L), eps);

Nt = numel(Trange)-1;
u = u0;
ll = zeros(Nt, 3);
tic;
for t=1:Nt
    d0t = d0ts(t);
    tau = (Trange(t+1)-Trange(t))/(2^d0t);
    
    Grad_t = IpaS(d0t,-1)/tau;
%     CN_term = IpaS(d0t,1)*0.5;   % Crank-Nicolson
    CN_term = tt_eye(2,d0t);   % Backw Euler
    e1t = tt_tensor(num2cell([1;0]*ones(1,d0t), 1));
    et = tt_ones(2, d0t);
    
    M = tkron(I, Grad_t)+tkron(A, CN_term);
    
    rhs = u/tau; %-0.5*A*u;
    rhs = tkron(rhs, e1t);
    U = tkron(u, et);
    
    U = amen_solve2(M, rhs, tol, 'x0', U);
    
    % Extract the final solution
    ext = tt_unit(2,d0t,2*ones(d0t,1));
    u = dot(ext, U, 3*L+1, U.d);
    
    mass = dot(u, tt_ones(u.n))
    u = u/mass;
    
    u2 = tt_reshape(u, 2^L*ones(1,3));
    % It should be the delta-function
    % Find the final distribution
    [v,l]=tt_stat(u2, 'lm');
    l=l-1
    ll(t,:)=l;
end
toc;

sp = zeros(2^L, 4);
for z=1:2^L
    [v,l]=tt_stat(u2(:,:,z), 'lm');
    if (v>tol)
        sp(z, 1)=l(1)-1;
        sp(z, 2)=l(2)-1;
        sp(z, 3)=z-1;
        sp(z, 4)=log2(v);
    else
        sp(z,:)=NaN;
    end
end
z = z-1;
scatter3(sp(1:z,1), sp(1:z,2), sp(1:z,3), 100, sp(1:z,4), 'filled');
xlabel('H');
ylabel('O');
zlabel('H2O');
title('log2(P)');
colorbar;
