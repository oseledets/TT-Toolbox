L = 7;

Trange = [0:5:20];
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

w1 = 2*mkron(o-e1-e2, o-e1, o);
w2 = 1*mkron(o, o, o-e1);
% w = mkron(x, x, o);

% CME matrix
A = -(mkron(S2,S1,S2')*diag(w1) - diag(w1)); % forward
A = A - (mkron(S2',S1',S2)*diag(w2) - diag(w2)); % backward
A = round(A, eps);

% initial guess
e1 = zeros(2^L,1); e1(l(1)+1)=1;
u0 = tt_tensor(e1);
e1 = zeros(2^L,1); e1(l(2)+1)=1;
u0 = kron(u0, tt_tensor(e1));
e1 = zeros(2^L,1); e1(l(3)+1)=1;
u0 = kron(u0, tt_tensor(e1));
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
%     CN_term = IpaS(d0t,1)*0.5;        
    CN_term = tt_eye(2,d0t);
    e1t = tt_tensor(num2cell([1;0]*ones(1,d0t), 1));
    et = tt_ones(2, d0t);
    
    M = kron(I, Grad_t)+kron(A, CN_term);
    
    rhs = u/tau; %-0.5*A*u;
    rhs = kron(rhs, e1t);
    U = kron(u, et);
    
    U = als_adapt_solve(M, rhs, tol, 'x0', U, 'max_full_size', 1, 'dirfilter', 1, 'local_prec', '', 'kicktype', 'tail', 'kickrank', 3, 'nswp', 15);
    
    % Extract the final solution
    ind = num2cell([1;2]*ones(1,3*L+d0t), 1);
    for i=3*L+1:3*L+d0t
        ind{i}=2;
    end;
    u = U(ind);
    u = tt_reshape(u, 2*ones(1,3*L));
    
    mass = dot(u, tt_ones(u.n))
    u = u/mass;
    
    u2 = tt_reshape(u, 2^L*ones(1,3));
    % It should be the delta-function
    % Find the final distribution
    [v,l]=tt_max_abs(u2);
    l=l-1
    ll(t,:)=l;
end;
toc;

sp = zeros(2^L, 4);
for z=1:2^L
    [v,l]=tt_max_abs(u2(:,:,z));
    if (v>tol)
        sp(z, 1)=l(1)-1;
        sp(z, 2)=l(2)-1;
        sp(z, 3)=z-1;
        sp(z, 4)=log2(v);
    else
        sp(z,:)=NaN;
    end;
end;
z = z-1;
% [v,l]=tt_max_abs(u2);
% z = l(3);
scatter3(sp(1:z,1), sp(1:z,2), sp(1:z,3), 100, sp(1:z,4), 'filled');

[V,L]=amr_eig_full(A, u0, tol, 5, 8);