L = 7;

Trange = [0:5e-3:0.1];
d0ts = 12*ones(1,numel(Trange)-1);

tol = 1e-5;
eps = 1e-12;

% initial distribution
l0 = [90,30,0]; % H2, N2, HH3
l=l0;

% One shift z=-1 (for N2)
S1 = tt_shf(L);
S2 = round(S1*S1, eps); %  double shift z=-2 (for NH3)
S3 = round(S2*S1, eps); %  triple shift z=-3 (for H2)
% "-z" corresponds to Sz'

I = tt_eye(2,3*L);

x = tt_x(2,L);
o = tt_ones(2,L);

damp_f1 = round(o-tt_unit(2,L,1)-tt_unit(2,L,2)-tt_unit(2,L,3), eps);
damp_f2 = round(o-tt_unit(2,L,1), eps);
damp_b = round(o-tt_unit(2,L,1)-tt_unit(2,L,2), eps);

w1 = 1.6e-4*mkron(damp_f1.*(x.*x.*x), damp_f2.*(x), o); % 400K
w2 = mkron(o, o, damp_b.*(x.*x));
% w = mkron(x, x, o);

% CME matrix
A = -(mkron(S3,S1,S2')*diag(w1) - diag(w1)); % forward
A = A - (mkron(S3',S1',S2)*diag(w2) - diag(w2)); % backward
A = round(A, eps);

% initial guess
u0 = mkron(tt_unit(2,L,l(1)+1), tt_unit(2,L,l(2)+1), tt_unit(2,L,l(3)+1));
% e1 = zeros(2^L,1); e1(l(1)+1)=1;
% u0 = tt_tensor(e1);
% e1 = zeros(2^L,1); e1(l(2)+1)=1;
% u0 = kron(u0, tt_tensor(e1));
% e1 = zeros(2^L,1); e1(l(3)+1)=1;
% u0 = kron(u0, tt_tensor(e1));
% % make QTT
% u0 = tt_reshape(u0, 2*ones(1,3*L), eps);

Nt = numel(Trange)-1;
meanconc = zeros(Nt, 3);
u = u0;
ll = zeros(Nt, 3);
tic;
for t=1:Nt
    d0t = d0ts(t);
    tau = (Trange(t+1)-Trange(t))/(2^d0t);
    
    Grad_t = IpaS(d0t,-1)/tau;
    CN_term = IpaS(d0t,1)*0.5;        
%     CN_term = tt_eye(2,d0t);
    e1t = tt_tensor(num2cell([1;0]*ones(1,d0t), 1));
    et = tt_ones(2, d0t);
    
    M = kron(I, Grad_t)+kron(A, CN_term);
    
    rhs = u/tau-0.5*A*u;
    rhs = kron(rhs, e1t);
    U = kron(u, et);
    
    U = als_adapt_solve(M, rhs, tol, 'x0', U, 'max_full_size', 1, 'dirfilter', 1, ...
        'local_prec', 'ljacobi', 'kicktype', 'one', 'kickrank', 3, ...
        'nswp', 10, 'resid_damp', 10, 'tol2', tol/2);
    
    % Extract the final solution
    ind = num2cell([1;2]*ones(1,3*L+d0t), 1);
    for i=3*L+1:3*L+d0t
        ind{i}=2;
    end;
    u = U(ind);
    u = tt_reshape(u, 2*ones(1,3*L));
    
    mass = dot(u, tt_ones(u.n))
    u = u/mass;
    
    meanconc(t,1) = dot(u, mkron(x,o,o));
    meanconc(t,2) = dot(u, mkron(o,x,o));
    meanconc(t,3) = dot(u, mkron(o,o,x));
    
    figure(1);
    plot(meanconc(1:t,:));
    legend('[H_2]', '[N_2]', '[NH_3]');
    drawnow;
    
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
        sp(z,:)=-Inf;
    end;
end;
z = z-1;
% [v,l]=tt_max_abs(u2);
% z = l(3);
figure(2)
scatter3(sp(1:z,1), sp(1:z,2), sp(1:z,3), 100, sp(1:z,4), 'filled');
colorbar;

outcome = meanconc(end, 3)/l0(2)*50

figure(3);
plot(2.^(sp(:,4)))

% [V,L]=amr_eig_full(A, u0, tol, 5, 8);