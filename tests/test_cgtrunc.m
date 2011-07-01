d = 7; n=2^d;
r = 9;
maxit_inner = 50;
maxit_outer = 100;

% % generate some tensor
% M = tt_qlaplace_dd([d,d]);
% M = tt_matrix(M)*(2^d+1)^2;
% f = tt_tensor(tt_ones(2*d,2));
% u = dmrg_solve2(M,f,f,1e-10,[], [], 10);
% resid = norm(M*u-f)/norm(f)
% 
% pause(1);
% 
% u = core(u);
% u = qtt_to_full(u, [n,n]);

% % add some shit
% u2 = u;
% u2(:,1:n/2)=rand(n,n/2);
% 
% 
% v = ones(n,n);
% v(:,1:n/2)=zeros(n,n/2);

[x,y]=meshgrid(-10+20*(1:1:n)/(n+1));
u2 = exp(-0.5*(x-y).^2);

v = exp(-0.5*(x.^2+y.^2));
W = full_to_tt(v, 1e-10);
W = tt_vec_to_diag(W);


[X,Y] = cg_trunc(W, u2, r, 1e-8, maxit_inner);
% We have to restart cg sometimes
for it=1:maxit_outer
    tic;
    [X,Y] = cg_trunc(W, u2, r, 1e-8, maxit_inner, X, Y);
    toc;
    fprintf('\n --restart--\n');
    pause(0.5);
end;