dx = 9;
dt = 15;
alpha = 1;

a = -10; b = 10;
T = 10;

h = (b-a)/(2^dx-1);
tau = T/(2^dt-1);

eps = 1e-6;

u0 = exp(-((a:h:b).^2)');
u0 = tt_tensor(full_to_qtt(u0, eps));
f = -exp(-(((a:h:b)+5).^2)')+0.5;
% f = tt_tensor(tt_ones(dx,2))*0;
f = tt_tensor(full_to_qtt(f, eps));

%%%%%%%%%% Spacial matrices
% Neumann lapl
lp = nd_to_full(tt_qlaplace_dd(dx));
% lp(1,1)=1;
lp(2^dx,2^dx)=1;
lp = tt_matrix(full_to_nd(lp, 1e-10, dx));
lp = lp*alpha/(h^2);

% Gradient
A = spdiags(ones(2^dx,1)*[-1,0,1], [-1 0 1], 2^dx, 2^dx);
% A(1,1)=-1;
A(2^dx,2^dx)=1;
A=A/h;
A = tt_matrix(full_to_nd(full(A), 1e-10, dx));

% Mass
M = spdiags(ones(2^dx,1)*[1/3,1/3,1/3], [-1 0 1], 2^dx, 2^dx);
M = tt_matrix(full_to_nd(full(M), 1e-10, dx));

% Quadratic 3-tensor
Q = matmat2quadr(A,M);
Q = round(Q, eps);

Ix = tt_matrix(tt_eye(2,dx));

%%%%%%%%%% Parts of the global CN matrix
CNterm = tt_matrix(IpaS(dt,1))*0.5;
e1t = cell(dt,1);
for i=1:dt
    e1t{i}=[1;0];
end;
e1t = tt_tensor(e1t);
CNterm = CNterm - diag(e1t*0.5);
CNterm = round(CNterm, 1e-10);
Grad_t = tt_matrix(IpaS(dt,-1))/tau;

Q_CN = matmat2quadr(CNterm, CNterm);
Q_CN = round(Q_CN, 1e-10);


% Global RHS
fgl = kron(u0, e1t)/tau + kron(f, tt_tensor(tt_ones(dt,2))-e1t);
fgl = round(fgl, eps);

% Global matrix and 3-array
Bgl = kron(Ix, Grad_t) + kron(lp, CNterm);
Bgl = round(Bgl, 1e-10);

Qgl = kron(Q, Q_CN);


t0 = tic;
U = dmrg_quadr_solve(Bgl, Qgl, fgl, [], eps, [], [], 20);
dmrg_time = toc(t0)
% fprintf('');