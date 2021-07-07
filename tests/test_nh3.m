% Chemical master equation for 3*H2 + N2 <-> 2*NH3

L = 7;  % QTT dimension in state variables

Trange = 0:5e-3:5e-2;
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

w1 = 1.6e-4*mtkron(damp_f1.*(x.*x.*x), damp_f2.*(x), o); % 400K
w2 = mtkron(o, o, damp_b.*(x.*x));
% w = mtkron(x, x, o);

% CME matrix
A = -(mtkron(S3,S1,S2')*diag(w1) - diag(w1)); % forward
A = A - (mtkron(S3',S1',S2)*diag(w2) - diag(w2)); % backward
A = round(A, eps);

% initial state
u0 = mtkron(tt_unit(2,L,l(1)+1), tt_unit(2,L,l(2)+1), tt_unit(2,L,l(3)+1));
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
    CN_term = IpaS(d0t,1)*0.5;   % Crank-Nicolson
%     CN_term = tt_eye(2,d0t);   % Backw Euler
    e1t = tt_tensor(num2cell([1;0]*ones(1,d0t), 1));
    et = tt_ones(2, d0t);
    
    M = tkron(I, Grad_t)+tkron(A, CN_term);
    
    rhs = u/tau-0.5*A*u;
    rhs = tkron(rhs, e1t);
    U = tkron(u, et);

    U = amen_solve2(M, rhs, tol, 'x0', U);
    
    % Extract the final solution
    ext = tt_unit(2,d0t,2*ones(d0t,1));
    u = dot(ext, U, 3*L+1, U.d);
    
    mass = dot(u, tt_ones(u.n))
    u = u/mass;
    
    meanconc(t,1) = dot(u, mtkron(x,o,o));
    meanconc(t,2) = dot(u, mtkron(o,x,o));
    meanconc(t,3) = dot(u, mtkron(o,o,x));
    
    figure(1);
    plot(Trange(2:t+1), meanconc(1:t,:));
    legend('[H_2]', '[N_2]', '[NH_3]');
    drawnow;
    
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
        sp(z,:)=-Inf;
    end
end
z = z-1;
figure(2)
scatter3(sp(1:z,1), sp(1:z,2), sp(1:z,3), 100, sp(1:z,4), 'filled');
xlabel('H2');
ylabel('N2');
zlabel('NH3');
title('log2(P)');
colorbar;

outcome = meanconc(end, 3)/l0(2)*50

figure(3);
plot(2.^(sp(:,4)))
xlabel('NH3');
title('P');
