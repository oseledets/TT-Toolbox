% Run file to produce the CME experiment from the paper
%   S. Dolgov, D. Savostyanov,
%   "Alternating minimal energy methods for linear systems in higher
%   dimensions"

try
    maxNumCompThreads(1);
catch
    % Just skip. Sometimes if you specify -singleCompThread in the command
    % line, MATLAB will fail at maxNumCompThreads with scary, so tell him
    % it's okay.
end;

n = 64*ones(20,1);
T = 10;
Lt = 12;
tol = 1e-6;
kickrank = 4;
nswp = 20;
trunc_norm = 'fro';
verb = 3;

tol_gmres = 1e-3; % I would not set anything lower

d = numel(n);

Ap = cell(d,1);
% Cascadic MPO -- production reactions
Ap{1} = zeros(1,n(1),n(1),3);
Ap{1}(1,:,:,1) = eye(n(1));
Ap{1}(1,:,:,2) = diag((0:n(1)-1)./(5+(0:n(1)-1)));
Ap{1}(1,:,:,3) = (diag(ones(n(1)-1,1),-1)-eye(n(1)))*0.7*diag([ones(n(1)-1,1);0]);
for i=2:d-1
    Ap{i} = zeros(3,n(i),n(i),3);
    Ap{i}(1,:,:,1)=eye(n(i));
    Ap{i}(3,:,:,3)=eye(n(i));
    Ap{i}(1,:,:,2)=diag((0:n(i)-1)./(5+(0:n(i)-1)));
    Ap{i}(2,:,:,3)=(diag(ones(n(i)-1,1),-1)-eye(n(i)))*diag([ones(n(i)-1,1);0]);
end;
Ap{d} = zeros(3,n(d),n(d));
Ap{d}(2,:,:) = (diag(ones(n(d)-1,1),-1)-eye(n(d)))*diag([ones(n(d)-1,1);0]);
Ap{d}(3,:,:) = eye(n(d));
% Laplace MPO -- destruction reactions
Ad = cell(d,1);
Ad{1} = zeros(1,n(1),n(1),2);
Ad{1}(1,:,:,1) = (diag(ones(n(1)-1,1),1)-eye(n(1)))*diag((0:n(1)-1))*0.07;
Ad{1}(1,:,:,2) = eye(n(1));
for i=2:d-1
    Ad{i} = zeros(2,n(i),n(i),2);
    Ad{i}(1,:,:,1)=eye(n(i));
    Ad{i}(2,:,:,2)=eye(n(i));
    Ad{i}(2,:,:,1)=(diag(ones(n(i)-1,1),1)-eye(n(i)))*diag((0:n(i)-1))*0.07;
end;
Ad{d} = zeros(2,n(d),n(d),1);
Ad{d}(1,:,:) = eye(n(d));
Ad{d}(2,:,:) = (diag(ones(n(d)-1,1),1)-eye(n(d)))*diag((0:n(d)-1))*0.07;

% QTT-tize each dimension separately, it is more stable due to smaller size
Apq = [];
Adq = [];
for i=1:d
    r1 = size(Ap{i},1);
    r2 = size(Ap{i},4);
    Aloc = cell2core(tt_matrix, Ap(i));
    Aloc = tt_reshape(Aloc, factor(n(i))'*[1,1], 1e-13, r1, r2);
    Apq = tkron(Apq, Aloc);
    
    r1 = size(Ad{i},1);
    r2 = size(Ad{i},4);
    Aloc = cell2core(tt_matrix, Ad(i));
    Aloc = tt_reshape(Aloc, factor(n(i))'*[1,1], 1e-13, r1, r2);
    Adq = tkron(Adq, Aloc);
end;

% Now the final CME operator
A = Apq+Adq;
I = tt_eye(A.n);


% Global time scheme
tau = T/2^Lt;
Gt = IpaS(Lt,-1)/tau;
iGt = tt_qtoepl(tkron(tt_ones(2,Lt), tt_tensor([0;1])), Lt)*tau;
Mt = tt_eye(2,Lt);
e1t = tt_unit(2,Lt,1);
et = tt_ones(2,Lt);
eNt = tt_unit(2,Lt,2^Lt);

% Global matrix
B = tkron(I,tt_eye(2,Lt)) - tkron(A,iGt*Mt);


% Initial state
u = [];
for i=1:d
    u = tkron(u, tt_unit(factor(n(i))', numel(factor(n(i))), 1)); %QTT
end;
U0 = tkron(u,et);

% RHS
f = tkron(u, et);

% Symmetrized matrix and RHS
B2 = B'*B;
B2 = round(B2, 1e-13);
f2 = B'*f;
f2 = round(f2, 1e-13);

% Use MEX library. It is safe, since in the last TT-Toolbox it switches to
% pure Matlab automatically, if you did not compile the MEXes.
ismex = true;

% Solvers
% Reference solution
tic;
U_ex = amen_solve2(B, f, tol*1e-3, 'x0', U0, 'nswp', nswp*3, 'kickrank', kickrank, 'trunc_norm', trunc_norm, 'ismex', true, 'verb', 1, 'max_full_size', 50, 'local_restart', 50, 'local_iters', 2, 'resid_damp', 2);
toc;
    
% ALS(SD)
tic;
[U_alstpz,td_alstpz] = alstpz_solve(B,f,tol, 'x0',U0, 'max_full_size', 50, 'kickrank', kickrank, 'nswp', nswp, 'kicktype', 'svd', 'symm', false, 'ismex', ismex, 'verb', verb);
toc;
% ALS(SD)+SYMM
tic;
[U_alstpz_s,td_alstpz_s] = alstpz_solve(B2,f2,tol, 'x0',U0, 'max_full_size', 50, 'kickrank', kickrank, 'nswp', nswp, 'kicktype', 'svd', 'symm', false, 'ismex', ismex, 'verb', verb);
toc;
% AMEN+SVD
tic;
[U_amen_svd,td_amen_svd] = amen_solve2(B, f, tol, 'x0', U0, 'nswp', nswp, 'kicktype', 'svd', 'kickrank', kickrank, 'trunc_norm', trunc_norm, 'ismex', ismex, 'verb', verb, 'max_full_size', 50, 'local_restart', 50, 'local_iters', 2, 'resid_damp', 2);
toc;
% AMEN+SVD+SYMM
tic;
[U_amen_svd_s,td_amen_svd_s] = amen_solve2(B2, f2, tol, 'symm', false, 'x0', U0, 'nswp', nswp, 'kicktype', 'svd', 'kickrank', kickrank, 'trunc_norm', trunc_norm, 'ismex', ismex, 'verb', verb, 'max_full_size', 50, 'local_restart', 50, 'local_iters', 2, 'resid_damp', 2);
toc;
% AMEN+ALS
tic;
[U_amen_als,td_amen_als] = amen_solve2(B, f, tol, 'x0', U0, 'nswp', nswp, 'kickrank', kickrank, 'trunc_norm', trunc_norm, 'ismex', ismex, 'verb', verb, 'max_full_size', 50, 'local_restart', 50, 'local_iters', 2, 'resid_damp', 2);
toc;
% AMEN+ALS+SYMM
tic;
[U_amen_als_s,td_amen_als_s] = amen_solve2(B2, f2, tol, 'symm', false, 'x0', U0, 'nswp', nswp, 'kickrank', kickrank, 'trunc_norm', trunc_norm, 'ismex', ismex, 'verb', verb, 'max_full_size', 50, 'local_restart', 50, 'local_iters', 2, 'resid_damp', 2);
toc;
% DMRG
tic;
[U_dmrg,td_dmrg] = dmrg_solve3(B, f, tol, 'x0', U0, 'nswp', nswp, 'kickrank', 0, 'trunc_norm', trunc_norm, 'ismex', ismex, 'verb', verb, 'max_full_size', 50, 'local_restart', 50, 'local_iters', 2, 'resid_damp', 2, 'step_drank', 0, 'step_dpow', 0, 'min_dpow', 0.5, 'dirfilter', 1);
toc;
% DMRG+SYMM
tic;
[U_dmrg_s,td_dmrg_s] = dmrg_solve3(B2, f2, tol, 'symm', false, 'x0', U0, 'nswp', nswp, 'kickrank', 0, 'trunc_norm', trunc_norm, 'ismex', ismex, 'verb', verb, 'max_full_size', 50, 'local_restart', 50, 'local_iters', 2, 'resid_damp', 2, 'step_drank', 0, 'step_dpow', 0, 'min_dpow', 0.5, 'dirfilter', 1);
toc;

% TT-GMRES
[U_gmres,td_gmres] = tt_gmres(core(B), core(f), tol_gmres, 10, 15, tol_gmres, tol_gmres, [], [], [], [], verb);

% This is a time-dep problem => check the last snapshots
if (~isempty(whos('U_ex')))
    u_ex = chunk(U_ex, 1, f.d-Lt)*dot(eNt, chunk(U_ex, f.d-Lt+1, f.d)).';
    u_ex = tt_reshape(u_ex, n);
end;
fprintf('Snapshot errors: \n');
if (~isempty(whos('U_alstpz')))
    u_alstpz = chunk(U_alstpz, 1, f.d-Lt)*dot(eNt, chunk(U_alstpz, f.d-Lt+1, f.d)).';
    u_alstpz = tt_reshape(u_alstpz, n);
    fprintf('alstpz:\t\t%3.5e\n', norm(u_alstpz-u_ex)/norm(u_ex));
end;
if (~isempty(whos('U_alstpz_s')))
    u_alstpz_s = chunk(U_alstpz_s, 1, f.d-Lt)*dot(eNt, chunk(U_alstpz_s, f.d-Lt+1, f.d)).';
    u_alstpz_s = tt_reshape(u_alstpz_s, n);
    fprintf('alstpz_s:\t\t%3.5e\n', norm(u_alstpz_s-u_ex)/norm(u_ex));
end;
if (~isempty(whos('U_amen_svd')))
    u_amen_svd = chunk(U_amen_svd, 1, f.d-Lt)*dot(eNt, chunk(U_amen_svd, f.d-Lt+1, f.d)).';
    u_amen_svd = tt_reshape(u_amen_svd, n);
    fprintf('amen_svd:\t\t%3.5e\n', norm(u_amen_svd-u_ex)/norm(u_ex));
end;
if (~isempty(whos('U_amen_svd_s')))
    u_amen_svd_s = chunk(U_amen_svd_s, 1, f.d-Lt)*dot(eNt, chunk(U_amen_svd_s, f.d-Lt+1, f.d)).';
    u_amen_svd_s = tt_reshape(u_amen_svd_s, n);
    fprintf('amen_svd_s:\t\t%3.5e\n', norm(u_amen_svd_s-u_ex)/norm(u_ex));
end;
if (~isempty(whos('U_amen_als')))
    u_amen_als = chunk(U_amen_als, 1, f.d-Lt)*dot(eNt, chunk(U_amen_als, f.d-Lt+1, f.d)).';
    u_amen_als = tt_reshape(u_amen_als, n);
    fprintf('amen_als:\t\t%3.5e\n', norm(u_amen_als-u_ex)/norm(u_ex));
end;
if (~isempty(whos('U_amen_als_s')))
    u_amen_als_s = chunk(U_amen_als_s, 1, f.d-Lt)*dot(eNt, chunk(U_amen_als_s, f.d-Lt+1, f.d)).';
    u_amen_als_s = tt_reshape(u_amen_als_s, n);
    fprintf('amen_als_s:\t\t%3.5e\n', norm(u_amen_als_s-u_ex)/norm(u_ex));
end;
if (~isempty(whos('U_dmrg')))
    u_dmrg = chunk(U_dmrg, 1, f.d-Lt)*dot(eNt, chunk(U_dmrg, f.d-Lt+1, f.d)).';
    u_dmrg = tt_reshape(u_dmrg, n);
    fprintf('dmrg:\t\t%3.5e\n', norm(u_dmrg-u_ex)/norm(u_ex));
end;
if (~isempty(whos('U_dmrg_s')))
    u_dmrg_s = chunk(U_dmrg_s, 1, f.d-Lt)*dot(eNt, chunk(U_dmrg_s, f.d-Lt+1, f.d)).';
    u_dmrg_s = tt_reshape(u_dmrg_s, n);
    fprintf('dmrg_s:\t\t%3.5e\n', norm(u_dmrg_s-u_ex)/norm(u_ex));
end;

% Compare with KSL
Nksl = 50;
% Set the ranks to the proper valus returned by AMEN
u_ksl = chunk(U_amen_als, 1, f.d-Lt)*dot(eNt, chunk(U_amen_als, f.d-Lt+1, f.d)).';
u_ksl = round(u_ksl, tol*0.1);
r_ksl = u_ksl.r; 
tic;
u_ksl = u+0*tt_rand(u.n,u.d,r_ksl,-1);
for j=1:Nksl
    u_ksl = tt_ksl_ml(u_ksl, A, tt_zeros(u_ksl.n), T/Nksl);
end;
ttimes_ksl=toc;
u_ksl = tt_reshape(u_ksl, n);
fprintf('ksl:\t\t%3.5e\n', norm(u_ksl-u_ex)/norm(u_ex));
fprintf('CPU Time of the KSL: %g\n', ttimes_ksl);


% % Result history processing
err_alstpz = zeros(nswp,1);
err_alstpz_s = zeros(nswp,1);
err_amen_svd = zeros(nswp,1);
err_amen_svd_s = zeros(nswp,1);
err_amen_als = zeros(nswp,1);
err_amen_als_s = zeros(nswp,1);
err_dmrg = zeros(nswp,1);
err_dmrg_s = zeros(nswp,1);
% Measure the errors
for i=1:nswp
    if (~isempty(whos('td_alstpz')))
        err_alstpz(i) = norm(td_alstpz{2}{end,i}-U_ex)/norm(U_ex);
    end;
    if (~isempty(whos('td_alstpz_s')))
        err_alstpz_s(i) = norm(td_alstpz_s{2}{end,i}-U_ex)/norm(U_ex);
    end;
    if (~isempty(whos('td_amen_svd')))
        err_amen_svd(i) = norm(td_amen_svd{2}{end,i}-U_ex)/norm(U_ex);
    end;
    if (~isempty(whos('td_amen_svd_s')))
        err_amen_svd_s(i) = norm(td_amen_svd_s{2}{end,i}-U_ex)/norm(U_ex);
    end;
    if (~isempty(whos('td_amen_als')))
        err_amen_als(i) = norm(td_amen_als{2}{end,i}-U_ex)/norm(U_ex);
    end;
    if (~isempty(whos('td_amen_als_s')))
        err_amen_als_s(i) = norm(td_amen_als_s{2}{end,i}-U_ex)/norm(U_ex);
    end;
    if (~isempty(whos('td_dmrg')))
        err_dmrg(i) = norm(td_dmrg{2}{end-1,i*2-1}-U_ex)/norm(U_ex);
    end;
    if (~isempty(whos('td_dmrg_s')))
        err_dmrg_s(i) = norm(td_dmrg_s{2}{end-1,i*2-1}-U_ex)/norm(U_ex);
    end;
end;

% Prepare the data in the TikZ-readable form
dats = [(1:nswp)', td_alstpz{1}(end,:)', err_alstpz, td_alstpz_s{1}(end,:)', err_alstpz_s, ...
    td_amen_svd{1}(end,:)', err_amen_svd, td_amen_svd_s{1}(end,:)', err_amen_svd_s, ...
    td_amen_als{1}(end,:)', err_amen_als, td_amen_als_s{1}(end,:)', err_amen_als_s, ...
    td_dmrg{1}(end-1,1:2:end)', err_dmrg, td_dmrg_s{1}(end-1,1:2:end)', err_dmrg_s];
% dats layout:
%   iter(1)     t_alstpz(2)     e_alstpz(3)     t_alstpz_s(4)       e_alstpz_s(5)
%   t_as(6)     e_as(7)         t_as_s(8)       e_as_s(9)
%   t_aa(10)    e_aa(11)        t_aa_s(12)      e_aa_s(13)
%   t_d(14)     e_d(15)         t_d_s(16)       e_d_s(17)

% Process the GMRES output. It is different...
if (~isempty(whos('td_gmres')))
    iter_gmres = min([find(td_gmres{1}==0, 1), numel(td_gmres{1})+1])-1;
    err_gmres = zeros(iter_gmres, 1);
    for i=1:iter_gmres
        err_gmres(i) = norm(tt_tensor(td_gmres{2}{i})-U_ex)/norm(U_ex);
    end;
    dat_gmres = [(1:iter_gmres)', td_gmres{1}(1:iter_gmres)', err_gmres];
end;

% Draw 'em for humans
% iter
figure(1);
semilogy(dats(:,1), dats(:,[3,5,7,9,11,13,15,17]), dat_gmres(1:min(iter_gmres,nswp),1), dat_gmres(1:min(iter_gmres,nswp),3));
legend('alstz', 'alstz-s', 'amen-svd', 'amen-svd-s', 'amen-als', 'amen-als-s', 'dmrg', 'dmrg-s', 'gmres');
% time
figure(2);
loglog(dats(:,2),dats(:,3), dats(:,4),dats(:,5), ...
    dats(:,6),dats(:,7), dats(:,8),dats(:,9), ...
    dats(:,10),dats(:,11), dats(:,12),dats(:,13), ...
    dats(:,14),dats(:,15), dats(:,16),dats(:,17), ...
    dat_gmres(:,2), dat_gmres(:,3));

% % Uncomment this if you want to draw the data elsewhere
% save('cme20_err.dat', '-ascii', 'dats');
% save('cme20_err_gmres.dat', '-ascii', 'dat_gmres');

