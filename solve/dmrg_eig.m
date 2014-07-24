% DMRG eigenvalue solver.
%   function [x,theta,testdata]=dmrg_eig(A, tol, varargin)
% Tries to solve the lowest eigenvalue (ground state) problem A*x=theta*x
% using the DMRG iteration.
% With default parameters, reconstructs the genuine 2-site DMRG for one
% targeted state, without enrichments.
% However, it can also work as the "block 1-site DMRG" from ref. [*], where
% at least 2 states must be targeted, but only one block at a time can be
% considered. See parameters 'b' and 'numblocks' below.
%
% A must be a Hermitian matrix in the tt_matrix format,
% tol is the relative truncation/stopping threshold,
% varargin may contain a sequence of tuning parameters of the form 
% 'parameter1_name', parameter1_value, 'parameter2_name', parameter2_value
% and so on. Available parameters and default values are the following.
%   o b: Number of eigenpairs targeted (default: 1)
%   o numblocks: Number of TT blocks included into the current supercore.
%       Can be either 2 for 2-site dmrg, or 1 for 1-site dmrg (default: 2)
%   o nswp: maximal number of DMRG sweeps (default: 20)
%   o max_full_size: if the local system size (rx(i)*n(i)*rx(i+1)) is less
%       than max_full_size, the local problem is solved via eig, otherwise via
%       lobpcg (iteratively) (default: 50)
%   o local_iters: maximal number of the lobpcg iterations (default: 100)
%   o usemex: whether shall we use the optimized MEX library eig3d_primme
%       instead of lobpcg. It safely switches off automatically if
%       eig3d_primme is not found in the MATLAB path, as well as on complex
%       data. To obtain eig3d_primme, you need to compile it in the
%       TT-Toolbox/fmex directory, please follow instructions there (default: true)
%   o resid_damp: residual gap for the local solvers, i.e. the real
%       threshold for e.g. lobpcg is tol/(sqrt(d)*resid_damp) (default: 2)
%   o trunc_norm: truncation strategy: either the standard Frobenius norm
%       ('fro'), or the residual ('resid') (default: 'fro')
%   o rmax: maximal TT rank (bond dimension) allowed (default: Inf)
%   o tol_exit: stopping tolerance (default: tol)
%   o verb: verbosity level: silent (0), sweep info (1), block info (2) or
%       the full debug returned in the testdata output (3) (default: 1)
%   o kickrank: size of the random enrichment (default: 0)
%   o x0: Initial guess in the tt_tensor format (default: random rank-2)
%
% Output data are the solution x in the tt_tensor format, 
% the eigenvalue theta, and 
% (optionally) the history for debug (or convergence tests) purposes testdata.
%
%
% Implementation: Sergey Dolgov (sergey.v.dolgov@gmail.com)
% 
% [*] Dolgov, Khoromskij, Oseledets, Savostyanov,
% "Computation of extreme eigenvalues in higher dimensions using block
% tensor train format", Comp. Phys. Comm. 2014, http://dx.doi.org/10.1016/j.cpc.2013.12.017
%


function [x,theta,testdata]=dmrg_eig(A, tol, varargin)
% Number of eigenstates
b = 1;
% Threshold on the local size between direct-iterative local solving
max_full_size=50;
% Number of local iters
local_iters = 100;
% Whether shall we use the local solver from the MEX library eig3d_primme
usemex = true;
% Number of sweeps
nswp=20;
% Residual gap for the local solver
resid_damp = 2;
% Truncation strategy
%trunc_norm = 'resid';
trunc_norm = 'fro';
% Max rank
rmax = Inf;
% Exit tol
tol_exit = tol;
% Verb
verb=1;
% Random enrichment rank
kickrank = 0;
% Number of blocks to consider: 2 or 1
numblocks = 2;
% Initial guess
x=[];

for i=1:2:length(varargin)-1
    switch lower(varargin{i})
        case 'nswp'
            nswp=varargin{i+1};
        case 'rmax'
            rmax=varargin{i+1};
        case 'x0'
            x=varargin{i+1};
        case 'verb'
            verb=varargin{i+1};
        case 'resid_damp'
            resid_damp=varargin{i+1};
        case 'trunc_norm'
            trunc_norm=varargin{i+1};
        case 'kickrank'
            kickrank=varargin{i+1};
        case  'max_full_size'
            max_full_size=varargin{i+1};
        case  'local_iters'
            local_iters=varargin{i+1};
        case 'usemex'
            usemex=varargin{i+1};
        case 'tol_exit'
            tol_exit=varargin{i+1};
        case 'b'
            b = varargin{i+1};
        case 'numblocks'
            numblocks = varargin{i+1};            
        otherwise
            error('Unknown tuning parameter "%s"', varargin{i});
    end;
end;

% Extract the dimensions
if (~isa(A, 'tt_matrix'))
    error('A must be given in a tt_matrix class');
end;
d = A.d;
n = A.n;
% Initialize from random initial state if not passed otherwise
if (isempty(x))
    x = tt_rand(n, d, b, -1);
else
    if (~isa(x, 'tt_tensor'))
        error('x0 must be given in a tt_tensor class');
    end;
    if (x.d~=d)
        error('Inconsistent dim(A) and dim(x0)');
    end;
    if (any(x.n~=n))
        error('Inconsistent size(A) and size(x0)');
    end;
end;

% Disable MEX if it does not exist
if (usemex)&&(exist('eig3d_primme', 'file')<2)
    warning('MEX local solver is not found, disabled');
    usemex = false;
end;

% More housekeeping
ra = A.r;
rx = x.r;
if (rx(d+1)>1)
    % Leave only the first component from the initial guess
    x = x*eye(rx(d+1),1);
    rx(d+1) = 1;
end;
crA = core2cell(A);
crx = core2cell(x);

% Threshold for local problems
real_tol = (tol/sqrt(d))/resid_damp;

% Interfaces
phixax = cell(d+1,1); phixax{1}=1; phixax{d+1}=1;

% This is some convergence output for test purposes
testdata = cell(3,1);
testdata{1} = zeros(d-numblocks+1, nswp); % times
testdata{2} = zeros(d-numblocks+1, nswp, b); % evs
testdata{3} = zeros(d-numblocks+1, nswp); % local res or dx
t_dmrg_eig = tic;

% Presetup: compute the initial projections
for i=d:-1:2
    % QR (Gauge conditions) right-to-left
    ux = reshape(crx{i}, rx(i), n(i)*rx(i+1));
    [ux,vx]=qr(ux.', 0);
    cr2 = reshape(crx{i-1}, rx(i-1)*n(i-1), rx(i));
    cr2 = cr2*vx.';
    rx(i) = size(ux,2);
    crx{i} = reshape(ux.', rx(i), n(i), rx(i+1));
    crx{i-1} = reshape(cr2, rx(i-1), n(i-1), rx(i));
    
    % Pass the reduction to the left
    phixax(i) = rightreduce_matrix(phixax(i+1), crx{i}, crA(i), crx{i}, rx(i),n(i),rx(i+1), 1,ra(i),ra(i+1), rx(i),n(i),rx(i+1));
    if (usemex)&&(~isreal(phixax{i}))
        warning('Complex data detected, turning MEX local solver off');
        usemex = false;
    end;
end;
% Initial guess for EVs
theta = rightreduce_matrix(phixax(2), crx{1}, crA(1), crx{1}, rx(1),n(1),rx(1+1), 1,ra(1),ra(1+1), rx(1),n(1),rx(1+1));
theta = squeeze(theta{1});
if (usemex)&&(~isreal(theta))
    warning('Complex data detected, turning MEX local solver off');
    usemex = false;
end;
theta = diag(theta);
if (b>1)
    crx{1} = randn(rx(1), n(1), rx(2), b);
end;

% DMRG sweeps
swp = 1;
i = 1; % Start from the first block
dir = 1; % Current direction: forward(+1) or backward(-1)
% Meausure the errors and # of operations
max_dx = 0;
max_res = 0;
max_matvecs = 0;

while (swp<=nswp)
    % Extract the matrix parts, accelerate a plenty of iterations with them
    % Phi1: 1, rx'1, rx1, ra1
    % Phi2: ra2, rx'2, rx2, 1
    Phi1 = phixax(i); 
    A1 = crA(i);
    if (numblocks==1)
        % one-block DMRG, everything is ready
        Phi2 = phixax(i+1);
        rx1 = rx(i);
        nloc = n(i);
        rx2 = rx(i+1);
        ra1 = ra(i);
        ra2 = ra(i+1);
        % sol_prev. It is also an initial guess.
        sol_prev = reshape(crx{i}, rx(i)*n(i)*rx(i+1), b);
    else
        % Two-block DMRG, we have to merge some of two blocks into one.
        Phi2 = phixax(i+2);
        A2 = crA(i+1);
        if ((rx(i)*n(i))<=(n(i)*n(i+1)))&&((rx(i)*n(i))<=(n(i+1)*rx(i+2)))
            % Merge Phi1 and A1
            Phi1 = reshape(Phi1{1}, rx(i)*rx(i), ra(i));
            A1 = reshape(A1{1}, ra(i), n(i)*n(i)*ra(i+1));
            Phi1 = Phi1*A1;
            Phi1 = reshape(Phi1, rx(i), rx(i), n(i), n(i)*ra(i+1));
            Phi1 = permute(Phi1, [1,3,2,4]);
            Phi1 = {reshape(Phi1, rx(i)*n(i), rx(i)*n(i), ra(i+1))};
            A1 = A2;
            rx1 = rx(i)*n(i);
            nloc = n(i+1);
            rx2 = rx(i+2);
            ra1 = ra(i+1);
            ra2 = ra(i+2);
        elseif ((n(i)*n(i+1)<=(rx(i)*n(i))))&&((n(i)*n(i+1))<=(n(i+1)*rx(i+2)))
            % Merge A1 and A2
            A1 = reshape(A1{1}, ra(i)*n(i)*n(i), ra(i+1));
            A2 = reshape(A2{1}, ra(i+1), n(i+1)*n(i+1)*ra(i+2));
            A1 = A1*A2;
            A1 = reshape(A1, ra(i)*n(i), n(i), n(i+1), n(i+1)*ra(i+2));
            A1 = permute(A1, [1,3,2,4]);
            A1 = {reshape(A1, ra(i), n(i)*n(i+1), n(i)*n(i+1), ra(i+2))};
            rx1 = rx(i);
            nloc = n(i)*n(i+1);
            rx2 = rx(i+2);
            ra1 = ra(i);
            ra2 = ra(i+2);
        else
            % Merge A2 and Phi2
            Phi2 = reshape(Phi2{1}, ra(i+2), rx(i+2)*rx(i+2));
            A2 = reshape(A2{1}, ra(i+1)*n(i+1)*n(i+1), ra(i+2));
            Phi2 = A2*Phi2;
            Phi2 = reshape(Phi2, ra(i+1)*n(i+1), n(i+1), rx(i+2), rx(i+2));
            Phi2 = permute(Phi2, [1,3,2,4]);
            Phi2 = {reshape(Phi2, ra(i+1), n(i+1)*rx(i+2), n(i+1)*rx(i+2))};
            rx1 = rx(i);
            nloc = n(i);
            rx2 = n(i+1)*rx(i+2);
            ra1 = ra(i);
            ra2 = ra(i+1);
        end;
        % sol_prev. It is also an initial guess.
        if (dir>0)
            sol_prev = reshape(crx{i}, rx(i)*n(i)*rx(i+1), b);
            sol_prev = sol_prev.';
            sol_prev = reshape(sol_prev, b*rx(i)*n(i), rx(i+1));
            sol_prev = sol_prev*reshape(crx{i+1}, rx(i+1), n(i+1)*rx(i+2));
            sol_prev = reshape(sol_prev, b, rx(i)*n(i)*n(i+1)*rx(i+2));
            sol_prev = sol_prev.';
        else
            sol_prev = reshape(crx{i+1}, rx(i+1), n(i+1)*rx(i+2)*b);
            sol_prev = reshape(crx{i}, rx(i)*n(i), rx(i+1))*sol_prev;
            sol_prev = reshape(sol_prev, rx(i)*n(i)*n(i+1)*rx(i+2), b);
        end;
    end;
    
    % Initial residual
    res_prev = norm(local_matvec(sol_prev, rx1,nloc,rx2,b, rx1,nloc,rx2, Phi1, A1, Phi2, 1,ra1,ra2)-sol_prev*diag(theta));
    
    % Check whether we actually need to solve anything
    if (res_prev>real_tol)
        if (rx1*nloc*rx2<max_full_size) % Full solution
            %      |     |    |
            % B = Phi1 - A1 - Phi2
            %      |     |    |
            Bxx = assemble_local_matrix(Phi1, A1, Phi2, 1,ra1,ra2, rx1,nloc,rx2, rx1,nloc,rx2);
            Bxx = (Bxx+Bxx')*0.5; % Ensure the symmetry. At least now...
            [sol,L]=eig(Bxx);
            L = diag(L);
            L = real(L); % exclude possible i*1e-15 noise
            [~,ind] = sort(L, 'ascend');
            theta = L(ind(1:b));
            sol = sol(:,ind(1:b));
            num_matvecs = 1;
        else % Iterative solution
            if (usemex)
                Phi1m = reshape(Phi1{1}, rx1, rx1, ra1);
                Phi2m = reshape(Phi2{1}, ra2*rx2, rx2);
                Phi2m = Phi2m.';
                Phi2m = reshape(Phi2m, rx2, ra2, rx2);
                A1m = reshape(A1{1}, ra1, nloc, nloc, ra2);
                [sol,theta, num_matvecs]=eig3d_primme(Phi1m,A1m,Phi2m,real_tol,b,sol_prev,local_iters);
            else
                [sol,theta,~,lambdaHistory]=lobpcg(sol_prev,@(x)local_matvec(x, rx1,nloc,rx2,[], rx1,nloc,rx2, Phi1, A1, Phi2, 1,ra1,ra2),  real_tol, local_iters);
                num_matvecs = numel(lambdaHistory);
            end;
        end;
        
        % count the number of MatVecs
        max_matvecs = max(max_matvecs, num_matvecs);
        
        if (~strcmp(trunc_norm, 'fro'))
            % We need the new residual for the corresp. trunc. strategy
            res_new = norm(local_matvec(sol, rx1,nloc,rx2,b, rx1,nloc,rx2, Phi1, A1, Phi2, 1,ra1,ra2)-sol*diag(theta));
        end;
    else
        res_new = res_prev;
        sol = sol_prev;
    end;

    % L2-norm convergence check. 
    dx = norm(sol*(sol'*sol_prev)-sol_prev);
    max_dx = max(max_dx, dx);
    max_res = max(max_res, res_prev);
    
    if ((dir>0)&&(i<d))||((dir<0)&&(i>1))||(numblocks==2) % SVD and enrichment
        if (dir>0)
            if (numblocks==2)
                % In 2-block DMRG, we have to separate modes
                sol = reshape(sol, rx(i)*n(i), n(i+1)*rx(i+2)*b);
            else
                sol = reshape(sol, rx(i)*n(i), rx(i+1)*b);
            end;
        else
            sol = sol.';
            if (numblocks==2)
                % In 2-block DMRG, we have to separate modes
                sol = reshape(sol, b*rx(i)*n(i), n(i+1)*rx(i+2));
            else
                sol = reshape(sol, b*rx(i), n(i)*rx(i+1));
            end;
        end;
        
        [ux,s,vx]=svd(sol, 'econ');
        s = diag(s);
        
        if (strcmp(trunc_norm, 'fro')) % Fro-norm truncation
            r = my_chop2(s, real_tol*resid_damp*norm(s));
        else
            % Residual truncation
            % start from the old rank
            r = min([rx(i),rx(i+1),numel(s)]);
            cursol = ux(:,1:r)*diag(s(1:r))*vx(:,1:r)';
            if (dir>0)
                cursol = reshape(cursol, rx1*nloc*rx2, b);
            else
                cursol = reshape(cursol, b, rx1*nloc*rx2);
                cursol = cursol.';
            end;
            res = norm(local_matvec(cursol, rx1,nloc,rx2,b, rx1,nloc,rx2, Phi1, A1, Phi2, 1,ra1,ra2)-cursol*diag(theta));
            if (res<max(real_tol, res_new)*resid_damp)
                drank = -1; % rank is overestimated; decrease
            else
                drank = 1; % residual is large; increase the rank
            end;
            while (r>0)&&(r<=numel(s))
                cursol = ux(:,1:r)*diag(s(1:r))*vx(:,1:r)';
                if (dir>0)
                    cursol = reshape(cursol, rx1*nloc*rx2, b);
                else
                    cursol = reshape(cursol, b, rx1*nloc*rx2);
                    cursol = cursol.';
                end;
                res = norm(local_matvec(cursol, rx1,nloc,rx2,b, rx1,nloc,rx2, Phi1, A1, Phi2, 1,ra1,ra2)-cursol*diag(theta));
                if (drank>0)
                    if (res<max(real_tol, res_new)*resid_damp)
                        break;
                    end;
                else
                    if (res>=max(real_tol, res_new)*resid_damp)
                        break;
                    end;
                end;
                r = r+drank;
            end;
            if (drank<0)
                r=r+1;
            end;
        end;
        
        r = min(r, numel(s));
        r = min(r, rmax);
        
        if (verb==2)
            fprintf('=dmrg_eig= swp=%d, block=%d, dx=%3.3e, r=%d\n', swp, i, dx, r);
        end;
    end;
    
    if (dir>0)&&(i<d) % Forward sweep
        ux = ux(:,1:r);
        vx = conj(vx(:,1:r))*diag(s(1:r));
        
        if (kickrank>0) 
            % Enrich the solution
            zx = randn(rx(i)*n(i), kickrank);
            % Concatenate the bases and reinforce the gauges.
            % In future: insert symmetries here?
            [ux,rv]=qr([ux,zx], 0);
            rv = rv(:,1:r);
            vx = vx*rv.';
        end;
        
        if (numblocks==1)
            % Drop the [r x r] factor to the next core -- White's prediction
            r = size(ux, 2);
            cr2 = crx{i+1};
            cr2 = reshape(cr2, rx(i+1), n(i+1)*rx(i+2));
            vx = reshape(vx, rx(i+1), b*r);
            cr2 = vx.'*cr2;
            cr2 = reshape(cr2, b, r*n(i+1)*rx(i+2));
            cr2 = cr2.';
        else
            % Replace the next core by vx
            cr2 = vx.';            
            cr2 = reshape(cr2, r, n(i+1), rx(i+2), b);
        end;
        
        % Stuff them back
        rx(i+1) = r;
        crx{i} = reshape(ux, rx(i), n(i), rx(i+1));
        crx{i+1} = reshape(cr2, rx(i+1), n(i+1), rx(i+2), b);
        
        % Compute new reductions
        phixax(i+1) = leftreduce_matrix(phixax(i), ux, crA(i), ux, rx(i),n(i),rx(i+1), 1,ra(i),ra(i+1), rx(i),n(i),rx(i+1));       
    elseif (dir<0)&&((i>1)||(numblocks==2)) % Backward sweep
        vx = conj(vx(:,1:r));
        ux = ux(:,1:r)*diag(s(1:r));
        
        if (kickrank>0)
            % Enrich the solution
            if (numblocks==1)
                zx = randn(n(i)*rx(i+1), kickrank);
            else
                zx = randn(n(i+1)*rx(i+2), kickrank);
            end;
            % Concatenate the bases and reinforce the gauges.
            % In future: insert symmetries here?
            [vx,rv]=qr([vx,zx], 0);
            rv = rv(:,1:r);
            ux = ux*rv.';
        end;
        
        if (numblocks==1)
            % Drop the [r x r] factor to the next core -- White's prediction
            r = size(vx, 2);
            cr2 = crx{i-1};
            cr2 = reshape(cr2, rx(i-1)*n(i-1), rx(i));
            ux = reshape(ux, b, rx(i)*r);
            ux = ux.';
            ux = reshape(ux, rx(i), r*b);
            cr2 = cr2*ux;
            cr2 = reshape(cr2, rx(i-1), n(i-1), r, b);            
            % Stuff them back
            rx(i) = r;
            vx = vx.';
            crx{i} = reshape(vx, rx(i), n(i), rx(i+1));
            crx{i-1} = reshape(cr2, rx(i-1), n(i-1), rx(i), b);            
            % Compute new reductions
            phixax(i) = rightreduce_matrix(phixax(i+1), vx, A1, vx, rx(i),n(i),rx(i+1), 1,ra(i),ra(i+1), rx(i),n(i),rx(i+1));            
        else
            % Replace both cores
            vx = vx.';
            vx = reshape(vx, r, n(i+1), rx(i+2));
            ux = reshape(ux, b, rx(i)*n(i)*r);
            ux = ux.';            
            rx(i+1) = r;
            crx{i+1} = vx;
            crx{i} = reshape(ux, rx(i), n(i), r, b);
            % Compute new reductions
            phixax(i+1) = rightreduce_matrix(phixax(i+2), vx, crA(i+1), vx, rx(i+1),n(i+1),rx(i+2), 1,ra(i+1),ra(i+2), rx(i+1),n(i+1),rx(i+2));
        end;
    else
        % Just stuff back the last core
        sol = reshape(sol, rx(i), n(i), rx(i+1), b);
        crx{i} = sol;
    end;
    
    if (verb>2)
        % Report the debug data if necessary
        if (dir>0)
            testdata{1}(i,swp) = toc(t_dmrg_eig);
            testdata{2}(i,swp,:) = theta;
            if (strcmp(trunc_norm, 'fro'))
                testdata{3}(i,swp) = dx;
            else
                testdata{3}(i,swp) = res_prev;
            end;
        else
            testdata{1}(d-numblocks+2-i,swp) = toc(t_dmrg_eig);
            testdata{2}(d-numblocks+2-i,swp,:) = theta;
            if (strcmp(trunc_norm, 'fro'))
                testdata{3}(d-numblocks+2-i,swp) = dx;
            else
                testdata{3}(d-numblocks+2-i,swp) = res_prev;
            end;
        end;
    end;
    
    i = i+dir;
    
    % Check for the end of the sweep
    if ((dir>0)&&(i>d-numblocks+1))||((dir<0)&&(i<1))
        % Report all
        if (verb>0)
            fprintf('=dmrg_eig= sweep %d, max_dx: %3.3e, max_res: %3.3e, erank: %g, theta: %3.15e, mv: %d\n', swp, max_dx, max_res, sqrt(rx(1:d)'*(n.*rx(2:d+1))/sum(n)), sum(theta), max_matvecs);
        end;
        
        % Check the stops
        if (strcmp(trunc_norm, 'fro'))
            if (max_dx<tol_exit)&&(verb<3)&&(dir>0)
                break;
            end;
        else
            if (max_res<tol_exit)&&(verb<3)&&(dir>0)
                break;
            end;
        end;
        
        swp = swp+1;
        % Go backward
        dir = -dir;
        i = i+dir;
        % Clear the errors
        max_dx = 0;
        max_res = 0;
        max_matvecs = 0;
    end;
end;

crx{d} = reshape(crx{d}, rx(d), n(d), b);
x = cell2core(tt_tensor, crx);

end


% Accumulates the left reduction W{1:k}'*A{1:k}*X{1:k}
function [WAX2] = leftreduce_matrix(WAX1, w, A, x, rw1,n,rw2, Ra,ra1,ra2, rx1,m,rx2)
% Left WAX has the form of the first matrix TT block, i.e. [rw, rx, ra]
WAX2 = WAX1;
wc = reshape(w, rw1, n*rw2);
xc = reshape(x, rx1*m, rx2);
for k=1:Ra
    WAX2{k} = reshape(WAX2{k}, rw1, rx1*ra1(k));
    WAX2{k} = wc'*WAX2{k}; % size n rw2 x rx1 ra1
    WAX2{k} = reshape(WAX2{k}, n, rw2*rx1*ra1(k));
    WAX2{k} = WAX2{k}.';
    WAX2{k} = reshape(WAX2{k}, rw2*rx1, ra1(k)*n);
    tmp = reshape(A{k}, ra1(k)*n, m*ra2(k));
    WAX2{k} = WAX2{k}*tmp; % size rw2 rx1 m ra2
    WAX2{k} = reshape(WAX2{k}, rw2, rx1*m*ra2(k));
    WAX2{k} = WAX2{k}.';
    WAX2{k} = reshape(WAX2{k}, rx1*m, ra2(k)*rw2);
    WAX2{k} = xc.'*WAX2{k}; % size rx2, ra2 rw2
    WAX2{k} = reshape(WAX2{k}, rx2*ra2(k), rw2);
    WAX2{k} = WAX2{k}.';
end;
end

% Accumulates the right reduction W{k:d}'*A{k:d}*X{k:d}
function [WAX1] = rightreduce_matrix(WAX2, w, A, x, rw1,n,rw2, Ra,ra1,ra2, rx1,m,rx2)
% Right WAX has the form of the last matrix TT block, i.e. [ra, rw, rx]
WAX1 = WAX2;
wc = reshape(w, rw1, n*rw2);
wc = conj(wc);
xc = reshape(x, rx1*m, rx2);
for k=1:Ra
    WAX1{k} = reshape(WAX1{k}, ra2(k)*rw2, rx2);
    WAX1{k} = xc*WAX1{k}.'; % size rx1 m x ra2 rw2
    WAX1{k} = reshape(WAX1{k}, rx1, m*ra2(k)*rw2);
    WAX1{k} = WAX1{k}.';
    WAX1{k} = reshape(WAX1{k}, m*ra2(k), rw2*rx1);
    tmp = reshape(A{k}, ra1(k)*n, m*ra2(k));
    WAX1{k} = tmp*WAX1{k}; % size ra1(k)*n, rw2*rx1
    WAX1{k} = reshape(WAX1{k}, ra1(k), n*rw2*rx1);
    WAX1{k} = WAX1{k}.';
    WAX1{k} = reshape(WAX1{k}, n*rw2, rx1*ra1(k));
    WAX1{k} = wc*WAX1{k}; % size rw1, rx1 ra1
    WAX1{k} = reshape(WAX1{k}, rw1*rx1, ra1(k));
    WAX1{k} = WAX1{k}.';
end;
end

% A matrix-vectors product for the matrix in the 3D TT (WAX1-A-WAX2), and
% full vectors of size (rx1*m*rx2) x b. Returns (rw1*n*rw2) x b
function [w]=local_matvec(x, rx1,m,rx2,b, rw1,n,rw2, WAX1, A, WAX2, Ra,ra1,ra2)
xc = reshape(x, rx1*m*rx2, []);
if (isempty(b))
    b = size(xc, 2);
end;
w = zeros(rw1*n*rw2, b);
xc = xc.';
xc = reshape(xc, b*rx1*m, rx2);
for k=1:Ra
    tmp = reshape(WAX2{k}, ra2(k)*rw2, rx2);
    wk = xc*tmp.';
    wk = reshape(wk, b*rx1, m*ra2(k)*rw2);
    wk = wk.';
    wk = reshape(wk, m*ra2(k), rw2*b*rx1);
    tmp = reshape(A{k}, ra1(k)*n, m*ra2(k));
    wk = tmp*wk;
    wk = reshape(wk, ra1(k)*n*rw2*b, rx1);
    wk = wk.';
    wk = reshape(wk, rx1*ra1(k), n*rw2*b);
    tmp = reshape(WAX1{k}, rw1, rx1*ra1(k));
    wk = tmp*wk;
    wk = reshape(wk, rw1*n*rw2, b);
    w = w+wk;
end;
end

% Builds the full (rw1*n*rw2) x (rx1*m*rx2) matrix from its TT blocks
function [B,sparseflag]=assemble_local_matrix(WAX1, A, WAX2, Ra,ra1,ra2, rw1,n,rw2, rx1,m,rx2)
% Check the sparsity of the matrix blocks
sparseflag = true;
for k=1:Ra
    if (~issparse(A{k}))
        sparseflag=false;
    end;
end;
if (sparseflag)
    B = sparse(rw2*rw1*n, rx2*rx1*m); % reverse order !!!
    % The reverse order is needed since usually the n x m part is large and
    % sparse, so let it be the senior dimension.
    % Note that currently only canonical sparse matrices are allowed
    for k=1:Ra
        tmp = reshape(WAX2{k}, rw2, rx2);
        tmp = sparse(tmp);
        Bk = reshape(WAX1{k}, rw1, rx1);
        Bk = sparse(Bk);
        Bk = kron(Bk, tmp); % mind endiannes
        Bk = kron(A{k}, Bk); % mind endiannes
        B = B+Bk;
    end;
else
    % There are dense blocks, everything is dense, and in the natural index
    % order
    B = zeros(rw1*n*rw2, rx1*m*rx2);
    for k=1:Ra
        Bk = reshape(WAX1{k}, rw1*rx1, ra1(k));
        tmp = reshape(A{k}, ra1(k), n*m*ra2(k));
        if (issparse(tmp))
            % Don't mix sparse if we are already full
            tmp = full(tmp);
        end;
        Bk = Bk*tmp;
        Bk = reshape(Bk, rw1, rx1, n, m*ra2(k));
        Bk = permute(Bk, [1,3,2,4]);
        Bk = reshape(Bk, rw1*n*rx1*m, ra2(k));
        tmp = reshape(WAX2{k}, ra2(k), rw2*rx2);
        Bk = Bk*tmp;
        Bk = reshape(Bk, rw1*n, rx1*m, rw2, rx2);
        Bk = permute(Bk, [1,3,2,4]);
        Bk = reshape(Bk, rw1*n*rw2, rx1*m*rx2);
        B = B+Bk;
    end;
end;
end

