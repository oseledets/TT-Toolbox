% Run file to produce the Poisson experiment from the paper
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


d = 16;
L = 6;

%!!! Last run: on dx360-11, E5-2670 0 @ 2.60GHz  !!!!!!!

tol = 1e-6;
nswp = 5;
kickrank = 4;

tol_gmres = 1e-3; % I would not set anything lower

% Matrix, RHS
A = tt_qlaplace_dd(L*ones(1,d))*(2^L+1)^2;
A = tt_reshape(A, 2^L*ones(d,2)); 
f = tt_ones(2^L, d);

% Reference solution
u_ex = amen_solve2(A, f, 1e-10, 'nswp', 30, 'max_full_size', 1500);

% Compare the solvers
% testdata format: {d,swp}
[u_svd,sd_svd] = amen_solve2(A, f, tol, 'x0', f, 'verb', 4, 'kicktype', 'svd', 'kickrank', kickrank, 'nswp', nswp, 'max_full_size', 50, 'trunc_norm', 'fro', 'ismex', true);
[u_als,sd_als] = amen_solve2(A, f, tol, 'x0', f, 'verb', 4, 'kicktype', 'als', 'kickrank', kickrank, 'nswp', nswp, 'max_full_size', 50, 'trunc_norm', 'fro', 'ismex', true);
[u_alstpz,sd_alstpz] = alstpz_solve(A,f,tol, 'x0',f, 'verb', 4, 'max_full_size', 50, 'kickrank', kickrank, 'nswp', nswp, 'kicktype', 'svd', 'symm', false, 'ismex', true);
[u_dmrg,sd_dmrg] = dmrg_solve3(A, f, tol, 'x0', f, 'kickrank', 0, 'min_dpow', 0.5, 'step_dpow', 0, 'step_drank', 0, 'max_full_size', 50, 'nswp', nswp, 'verb', 4, 'ismex', true, 'trunc_norm', 'fro', 'dirfilter', 1);
[U_gmres,td_gmres] = tt_gmres(core(A), core(f), tol_gmres, 4, 15, tol_gmres, tol_gmres, [], [], [], [], 3);


% Measure the A-errors
err0 = sqrt(dot(u_ex, A*u_ex));
errs_svd = zeros(d,nswp);
errs_als = zeros(d,nswp);
errs_alstpz = zeros(d,nswp);
errs_dmrg = zeros(d-1,nswp);
for i=1:nswp
    for j=1:d
        errs_svd(j,i) = sqrt(dot(sd_svd{2}{j,i}-u_ex, A*(sd_svd{2}{j,i}-u_ex)))/err0;
        errs_als(j,i) = sqrt(dot(sd_als{2}{j,i}-u_ex, A*(sd_als{2}{j,i}-u_ex)))/err0;
        errs_alstpz(j,i) = sqrt(dot(sd_alstpz{2}{j,i}-u_ex, A*(sd_alstpz{2}{j,i}-u_ex)))/err0;
        if (j<d)
            errs_dmrg(j,i) = sqrt(dot(sd_dmrg{2}{j,2*i-1}-u_ex, A*(sd_dmrg{2}{j,2*i-1}-u_ex)))/err0;
        end;
    end;
end;

% Prepare the data in the TikZ-readable form
%   1       2       3       4       5       6       7
% iter, t_alstpz, e_alstpz, t_svd, e_svd, t_als, e_als
dat_amens = [(1/d:1/d:nswp)', sd_alstpz{1}(:), errs_alstpz(:), sd_svd{1}(:), errs_svd(:), sd_als{1}(:), errs_als(:)];
dat_dmrg = [(1/(d-1):1/(d-1):nswp)', reshape(sd_dmrg{1}(1:d-1,1:2:nswp*2), [],1), errs_dmrg(:)];

% Process the GMRES output. It is different...
if (~isempty(whos('td_gmres')))
    iter_gmres = min([find(td_gmres{1}==0, 1), numel(td_gmres{1})+1])-1;
    err_gmres = zeros(iter_gmres, 1);
    for i=1:iter_gmres
        err_gmres(i) = sqrt(dot(tt_tensor(td_gmres{2}{i})-u_ex, A*(tt_tensor(td_gmres{2}{i})-u_ex)))/err0;
    end;
    dat_gmres = [(1:iter_gmres)', td_gmres{1}(1:iter_gmres)', err_gmres];
end;

% Draw 'em for humans
% iter
figure(1);
semilogy(dat_amens(:,1), dat_amens(:,[3,5,7]), dat_dmrg(:,1), dat_dmrg(:,3), dat_gmres(1:min(iter_gmres,nswp),1), dat_gmres(1:min(iter_gmres,nswp),3));
legend('alstpz', 'amen-svd', 'amen-als', 'dmrg', 'gmres');
% time
figure(2);
loglog(dat_amens(:,2), dat_amens(:,3), dat_amens(:,4), dat_amens(:,5), dat_amens(:,6), dat_amens(:,7), dat_dmrg(:,2), dat_dmrg(:,3), dat_gmres(:,2), dat_gmres(:,3));

% % Uncomment this if you want to draw the data elsewhere
% save('conv_lp16_amens.dat', '-ascii', 'dat_amens');
% save('conv_lp16_dmrg.dat', '-ascii', 'dat_dmrg');
% save('conv_lp16_gmres.dat', '-ascii', 'dat_gmres');
