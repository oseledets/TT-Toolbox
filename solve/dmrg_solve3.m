function [x,somedata]=dmrg_solve3(A, y, tol, varargin)
%Solution of linear systems in TT-format via DMRG iteration
%   [X,SWEEPS]=DMRG_SOLVE3(A,Y,TOL,OPTIONS) Attempts to solve the linear
%   system A*X = Y with accuracy/residual TOL using the two-sided DMRG iteration.
%   Matrix A has to be given in the TT-format, right-hand side Y should be
%   given in the TT-format also. Options are provided in form
%   'PropertyName1',PropertyValue1,'PropertyName2',PropertyValue2 and so
%   on. The parameters are set to default (in brackets in the following) 
%   The list of option names and default values are:
%     ------------------- Main options for DMRG in general
%        o nswp - maximal number of sweeps [20]
%        o rmax - maximal TT-rank of the solution [1000]
%        o x0 - initial guess [random rank-2 tensor]
%        o verb - verbosity level, 0-silent, 1-sweep info, 2-block info [1]
%        o dirfilter - shall we update the solution on the forward sweep
%          only (1), on backward sweep only (-1), or both (0) [0]
%        o symm - shall we symmetrize the problem (A'Ax=A'y) before 
%          the solution [false]
%        o trunc_norm - truncate in either Frob. ('fro'), or residual norm
%          ('residual') ['residual']
%        o ismex - shall we use the MEX lib solve3d_2 for local solution
%          instead of gmres. It safely switches off automatically if
%          solve3d_2 is not found in the MATLAB path, as well as on complex
%          data. To obtain solve3d_2, you need to compile it in the
%          TT-Toolbox/fmex directory, please follow instructions there [true]
%     ------------------- Options of the local solver
%        o local_restart - dimension of local gmres [40]
%        o local_iters - number of local gmres restarts [2]
%        o max_full_size - maximal size of the local matrix for the full solver [2500]
%        o resid_damp - solve local problems with accuracy tol/resid_damp.
%          Larger value may reduce a spurious noise from inexact local
%          solutions, but increase CPU time [2]
%        o local_prec - local preconditioner: '' (no prec.), 'ljacobi',
%          'cjacobi', 'rjacobi' ['']
%     ------------------- Heuristics that helped in some problems, but generally not clear
%        o kickrank - random enrichment size [2]
%        o min_dpow - sometimes it helps to truncate intermediate
%          iterations up to tol/(d^dpow) instead of tol/sqrt(d). This
%          parameter specifies the minimal dpow [1]
%        o step_dpow - a step of increasing/decreasing dpow, according to
%          the current convergence pattern [0.1]
%        o bot_conv - if the error changed by this factor or smaller, decrease
%          dpow by step_dpow and drank by step_drank [0.1]
%        o top_conv - if the error changed by this factor or larger,  increase
%          dpow by step_dpow and drank by step_drank [0.99]
%        o step_drank - a step of additional increasing/decreasing of ranks
%        o block_order - allows to update only a chunk of the solution. The
%          iteration starts from the first block, then goes for
%          block_order(1) blocks, then (from that position) for
%          block_order(2) blocks, and so on. For example, to update the
%          first k blocks, set block_order = [k,-k]. Normally this is not
%          needed, so the default value is [d-1, -(d-1)].
%
%       Example:
%           d=8; f=8; 
%           mat=tt_qlaplace_dd(d*ones(1,f)); %Laplace in the QTT-format
%           rhs=tt_ones(2,d*f); % Right-hand side of all ones
%           sol=dmrg_solve3(mat,rhs,1e-6);
%
%
% TT-Toolbox 2.2, 2009-2012
%
%This is TT Toolbox, written by Ivan Oseledets et al.
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------


% Inner parameters
max_full_size=2500;
step_dpow = 0.1; % stepsize for d-power in truncations
min_dpow = 1; % Minimal d-power for truncation
step_drank = 1;

resid_damp = 2; % Truncation error to true residual treshold
bot_conv = 0.1; % bottom convergence factor - if better, we can decrease dpow and drank
top_conv = 0.99; % top convergence factor - if worse, we have to increase dpow and drank


nswp=20;
local_restart=40;
local_iters=2;

local_prec = '';
% local_prec = 'jacobi';

rmax=1000;
trunc_norm = 'residual';
% trunc_norm = 'fro';

ismex = true;

% dirfilter = -1; % only backward
% dirfilter = 1; % only forward
dirfilter = 0; % both

symm = false;

verb=1;
kickrank = 2;
x=[];
block_order = [];

for i=1:2:length(varargin)-1
    switch lower(varargin{i})
        case 'nswp'
            nswp=varargin{i+1};
        case 'rmax'
            rmax=lower(varargin{i+1});
        case 'x0'
            x=varargin{i+1};
        case 'verb'
            verb=varargin{i+1};
        case 'local_prec'
            local_prec=varargin{i+1};
        case 'local_restart'
            local_restart=varargin{i+1};
        case 'local_iters'
            local_iters=varargin{i+1};
        case 'dirfilter'
            dirfilter=varargin{i+1};
        case 'kickrank'
            kickrank=varargin{i+1};
        case 'symm'
            symm=varargin{i+1};            
        case  'max_full_size'
            max_full_size=varargin{i+1};
        case 'step_dpow'
            step_dpow=varargin{i+1};
        case 'min_dpow'
            min_dpow=varargin{i+1};
        case 'step_drank'
            step_drank=varargin{i+1};            
        case 'resid_damp'
            resid_damp = varargin{i+1};
        case 'trunc_norm'
            trunc_norm = varargin{i+1};
        case 'bot_conv'
            bot_conv=varargin{i+1};
        case 'top_conv'
            top_conv=varargin{i+1};
        case 'block_order'
            block_order=varargin{i+1};    
        case 'ismex'
            ismex=varargin{i+1};            
            
        otherwise
            error('Unrecognized option: %s\n',varargin{i});
    end
end

local_prec_char = 0;
if (strcmp(local_prec, 'cjacobi')); local_prec_char = 1;  end;
if (strcmp(local_prec, 'ljacobi')); local_prec_char = 2;  end;
if (strcmp(local_prec, 'rjacobi')); local_prec_char = 3;  end;

if (A.n~=A.m)
    error(' DMRG does not know how to solve rectangular systems!\n Use dmrg_solve3(ctranspose(A)*A, ctranspose(A)*f, tol) instead.');
end;

% Disable MEX if it does not exist
if (ismex)&&(exist('solve3d_2', 'file')<2)
    warning('MEX local solver is not found, disabled');
    ismex = false;
end;

d = y.d;
n = A.n;
if (isempty(x))
    x = tt_rand(n, A.d, 2);
end;

if (isempty(block_order))
    block_order = [+(d-1), -(d-1)];
end;

if (symm)
    ry = (y.r).*(A.r);
    ra = (A.r).^2;
    crA = core2cell(A'*A);
    cry = core2cell(A'*y);
else
    ry = y.r;
    ra = A.r;
    crA = core2cell(A);
    cry = core2cell(y);
end;

rx = x.r;
crx = core2cell(x);

phia = cell(d+1,1); phia{1}=1; phia{d+1}=1;
phiy = cell(d+1,1); phiy{1}=1; phiy{d+1}=1;

somedata = cell(3,1);
somedata{1}=zeros(d,nswp*2);
somedata{2}=cell(d, nswp*2);
somedata{3}=zeros(d,nswp*2);

t_dmrg_solve = tic;

% Orthogonalization
for i=d:-1:2
    cr = crx{i};
    cr = reshape(cr, rx(i), n(i)*rx(i+1));
    [cr, rv]=qr(cr.', 0);
    cr2 = crx{i-1};
    cr2 = reshape(cr2, rx(i-1)*n(i-1), rx(i));
    cr2 = cr2*(rv.');
    rx(i) = size(cr, 2);
    cr = reshape(cr.', rx(i), n(i), rx(i+1));
    crx{i-1} = reshape(cr2, rx(i-1), n(i-1), rx(i));
    crx{i} = cr;
    
    phia{i} = compute_next_Phi(phia{i+1}, cr, crA{i}, cr, 'rl');
    phiy{i} = compute_next_Phi(phiy{i+1}, cr, [], cry{i}, 'rl');   
    
    if ((~isreal(phia{i}))||(~isreal(phiy{i})))&&(ismex)
        warning('Complex data detected, turning MEX local solver off');
        ismex = false;
    end;
end;


last_sweep = false;
swp = 1;
i = 1;

dx_old = ones(d,1);
dx = zeros(d,1);
max_res = 0;
max_dx = 0;
max_iter = 0;
% For extra-rank addition
dpows = ones(d,1)*min_dpow;
dranks = zeros(d,1);

cur_order = block_order;
order_index = 1;
dir = sign(cur_order(order_index));

% DMRG sweeps
while (swp<=nswp)
    % Extract elements - matrix
    Phi1 = phia{i}; Phi2 = phia{i+2};
    A1 = crA{i}; A2 = crA{i+1};
    % RHS
    rhs = phiy{i};
    rhs = rhs*reshape(cry{i}, ry(i), n(i)*ry(i+1));
    rhs = reshape(rhs, rx(i)*n(i), ry(i+1));
    rhs = rhs*reshape(cry{i+1}, ry(i+1), n(i+1)*ry(i+2));
    rhs = reshape(rhs, rx(i)*n(i)*n(i+1), ry(i+2));
    rhs = rhs*(phiy{i+2}.');
    rhs = reshape(rhs, rx(i)*n(i)*n(i+1)*rx(i+2),1);
    norm_rhs = norm(rhs);
    % sol_prev
    sol_prev = reshape(crx{i}, rx(i)*n(i), rx(i+1));
    sol_prev = sol_prev*reshape(crx{i+1}, rx(i+1), n(i+1)*rx(i+2));
    sol_prev = reshape(sol_prev, rx(i)*n(i)*n(i+1)*rx(i+2),1);
    
    real_tol = (tol/(d^dpows(i)))/resid_damp;
    if (last_sweep)
        real_tol = (tol/sqrt(d))/resid_damp;
    end;
    
    if (rx(i)*n(i)*n(i+1)*rx(i+2)<max_full_size) % Direct solution
        %      |     |    |     |
        % B = Phi1 - A1 - A2 - Phi2
        %      |     |    |     |
        if ((dir+dirfilter)~=0)
            B = reshape(permute(Phi1, [1, 3, 2]), rx(i)*rx(i), ra(i));
            B = B*reshape(A1, ra(i), n(i)*n(i)*ra(i+1));
            B = reshape(B, rx(i), rx(i), n(i), n(i), ra(i+1));
            B = permute(B, [1, 3, 2, 4, 5]);
            B = reshape(B, rx(i)*n(i)*rx(i)*n(i), ra(i+1));
            B = B*reshape(A2, ra(i+1), n(i+1)*n(i+1)*ra(i+2));
            B = reshape(B, rx(i)*n(i), rx(i)*n(i), n(i+1), n(i+1), ra(i+2));
            B = permute(B, [1, 3, 2, 4, 5]);
            B = reshape(B, rx(i)*n(i)*n(i+1)*rx(i)*n(i)*n(i+1), ra(i+2));
            B = B*reshape(permute(Phi2, [2, 1, 3]), ra(i+2), rx(i+2)*rx(i+2));
            B = reshape(B, rx(i)*n(i)*n(i+1), rx(i)*n(i)*n(i+1), rx(i+2), rx(i+2));
            B = permute(B, [1, 3, 2, 4]);
            B = reshape(B, rx(i)*n(i)*n(i+1)*rx(i+2), rx(i)*n(i)*n(i+1)*rx(i+2));
        end;
        
        res_prev = norm(bfun3(Phi1, A1, A2, Phi2, sol_prev) - rhs)/norm_rhs;

        if (res_prev>real_tol)&&((dir+dirfilter)~=0)
            sol = B \ rhs;
            flg = 0;
            res_new = norm(B*sol-rhs)/norm_rhs;
            iter=1;
        else
            sol = sol_prev;
            res_new = res_prev;
            flg=0;
            iter=0;
        end;
        
    else % Structured solution.
        
        res_prev = norm(bfun3(Phi1, A1, A2, Phi2, sol_prev) - rhs)/norm_rhs;
        
        if (res_prev>real_tol)&&((dir+dirfilter)~=0)
            % Since we only have procedures for 3 blocks, we have to
            % merge two of them. Usually this is not a serious
            % overhead... on those problems where DMRG works in
            % principle
            if (rx(i)*n(i)<=min(n(i)*n(i+1), n(i+1)*rx(i+2)))
                Phi1mex = reshape(permute(Phi1, [1, 3, 2]), rx(i)*rx(i), ra(i));
                Phi1mex = Phi1mex*reshape(A1, ra(i), n(i)*n(i)*ra(i+1));
                Phi1mex = reshape(Phi1mex, rx(i), rx(i), n(i), n(i), ra(i+1));
                Phi1mex = permute(Phi1mex, [1, 3, 2, 4, 5]);
                Phi1mex = reshape(Phi1mex, rx(i)*n(i), rx(i)*n(i), ra(i+1));
                Amex = reshape(A2, ra(i+1), n(i+1), n(i+1), ra(i+2));
                Phi2mex = permute(Phi2, [3,2,1]);
            elseif (n(i)*n(i+1)<=min(rx(i)*n(i), n(i+1)*rx(i+2)))
                Amex = reshape(A1, ra(i)*n(i)*n(i), ra(i+1));
                Amex = Amex*reshape(A2, ra(i+1), n(i+1)*n(i+1)*ra(i+2));
                Amex = reshape(Amex, ra(i), n(i), n(i), n(i+1), n(i+1), ra(i+2));
                Amex = permute(Amex, [1, 2, 4, 3, 5, 6]);
                Amex = reshape(Amex, ra(i), n(i)*n(i+1), n(i)*n(i+1), ra(i+2));
                Phi1mex = permute(Phi1, [1,3,2]);
                Phi2mex = permute(Phi2, [3,2,1]); % rx', ra, rx -> rx,ra,rx'
            else
                Phi1mex = permute(Phi1, [1,3,2]);
                Amex = reshape(A1, ra(i), n(i), n(i), ra(i+1));
                Phi2mex = reshape(A2, ra(i+1)*n(i+1)*n(i+1), ra(i+2));
                Phi2mex = Phi2mex*reshape(permute(Phi2, [2, 1, 3]), ra(i+2), rx(i+2)*rx(i+2));
                Phi2mex = reshape(Phi2mex, ra(i+1), n(i+1), n(i+1), rx(i+2), rx(i+2));
                Phi2mex = permute(Phi2mex, [3,5, 1, 2,4]);
                Phi2mex = reshape(Phi2mex, n(i+1)*rx(i+2), ra(i+1), n(i+1)*rx(i+2));
            end;
            
            if (~ismex)
                [sol,flg,iter] = solve3d_2ml(Phi1mex, Amex, Phi2mex, rhs, real_tol*norm_rhs, sol_prev, local_prec_char, local_restart, local_iters);
            else % ismex
                sol = solve3d_2(Phi1mex, Amex, Phi2mex, rhs, real_tol, 1, sol_prev, local_prec_char, local_restart, local_iters, 0);
                flg = 0;
                iter = 1;
            end;
            
            res_new = norm(bfun3(Phi1, A1, A2, Phi2, sol) - rhs)/norm_rhs;            
        else
            sol = sol_prev;
            res_new = res_prev;
            flg=0;
            iter=0;
        end;
    end;
    
    if (flg>0)
        fprintf('-warn- local solver did not converge at block %d\n', i);
    end;
    if (res_prev/res_new<resid_damp)&&(res_new>real_tol)&&((dir+dirfilter)~=0)
        fprintf('--warn-- the residual damp was smaller than in the truncation\n');
    end;
    
    dx(i) = norm(sol-sol_prev)/norm(sol);
    max_dx = max(max_dx, dx(i));
    
    if (strcmp(trunc_norm, 'fro'))
        somedata{3}(i,(swp-1)*2+1.5-dir/2) = dx(i);
    else
        somedata{3}(i,(swp-1)*2+1.5-dir/2) = res_prev;
    end;
    
    % The new core does not converge - increase rank
    if (dx(i)/dx_old(i)>top_conv)&&(dx(i)>tol)
        dranks(i)=dranks(i)+step_drank;
        dpows(i)=dpows(i)+step_dpow;
    end;
    % The new core converges well - try to decrease rank
    if (dx(i)/dx_old(i)<bot_conv)||(dx(i)<tol)
        dranks(i)=max(dranks(i)-step_drank, 0);
        dpows(i)=max(dpows(i)-step_dpow, min_dpow);
    end;
    
    if (last_sweep)
        dpows(i)=0.5;
        dranks(i)=0;
    end;
    
    % Error indicators
    max_res = max(max_res, res_prev);
    max_iter = max(max_iter, iter);
    
    % Truncation
    if ((dir+dirfilter)~=0)
        sol = reshape(sol, rx(i)*n(i), n(i+1)*rx(i+2));
        [u,s,v]=svd(sol, 'econ');
        s = diag(s);
        
        if (strcmp(trunc_norm, 'fro')) % We are happy with L2 truncation
            r = my_chop2(s, max(tol/(d^dpows(i)), res_new*resid_damp)*norm(s));
        else
            % Residual trunc; First, bin-search
            r1 = 1; r2 = numel(s); r = round((r1+r2)/2);
            while (r2-r1>1)
                cursol = u(:,1:r)*diag(s(1:r))*(v(:,1:r)');
                if (rx(i)*n(i)*n(i+1)*rx(i+2)<max_full_size)
                    res = norm(B*cursol(:)-rhs)/norm_rhs;
                else
                    res = norm(bfun3(Phi1, A1, A2, Phi2, cursol)-rhs)/norm_rhs;
                end;
                if (res<max(tol/(d^dpows(i)), res_new*resid_damp))
                    r2 = r;
                else
                    r1 = r;
                end;
                r = round((r1+r2)/2);
            end;
            % More accurate Linear search
            while (r<=numel(s))
                cursol = u(:,1:r)*diag(s(1:r))*(v(:,1:r)');
                if (rx(i)*n(i)*n(i+1)*rx(i+2)<max_full_size)
                    res = norm(B*cursol(:)-rhs)/norm_rhs;
                else
                    res = norm(bfun3(Phi1, A1, A2, Phi2, cursol)-rhs)/norm_rhs;
                end;
                if (res<max(tol/(d^dpows(i)), res_new*resid_damp))
                    break;
                end;
                r = r+1;
            end;
        end;
    else
        if (dir>0)
            [u,v]=qr(reshape(crx{i}, rx(i)*n(i), rx(i+1)), 0);
            v = v*reshape(crx{i+1}, rx(i+1), n(i+1)*rx(i+2));
            v = v';
            r = size(u,2);
            s = ones(r,1);
        else
            [v,u]=qr(reshape(crx{i+1}, rx(i+1), n(i+1)*rx(i+2)).', 0);
            v=conj(v);
            u=reshape(crx{i}, rx(i)*n(i), rx(i+1))*u.';
            r = size(u,2);
            s = ones(r,1);
        end;
    end;
    
    % Artificial rank increasing
    r = r+dranks(i);
    r = min(r, numel(s));
    r = min(r, rmax);
    
    if (verb==2)
        fprintf('=dmrg_solve3=   block %d{%d}, dx: %3.3e, iter: %d, r: %d\n', i, dir, dx(i), iter, r);
    end;
    
    if (dir>0) % left-to-right
        u = u(:,1:r);
        v = conj(v(:,1:r))*diag(s(1:r));
        % kick
        if (~last_sweep)&&((dir+dirfilter)~=0)
            [u,rv]=qr([u, randn(rx(i)*n(i), kickrank)], 0);
            v = [v, zeros(n(i+1)*rx(i+2), kickrank)];
            v = v*(rv.');
        end;
        r = size(u,2);
        u = reshape(u, rx(i), n(i), r);
        v = reshape(v.', r, n(i+1), rx(i+2));
        
        % Recompute phi. Left ones, so permute them appropriately
        phia{i+1} = compute_next_Phi(phia{i}, u, crA{i}, u, 'lr');
        phiy{i+1} = compute_next_Phi(phiy{i}, u, [], cry{i}, 'lr');
        
    else % right-to-left
        u = u(:,1:r)*diag(s(1:r));
        v = conj(v(:,1:r));
        % kick
        if (~last_sweep)&&((dir+dirfilter)~=0)
            [v,rv] = qr([v, randn(n(i+1)*rx(i+2), kickrank)], 0);
            u = [u, zeros(rx(i)*n(i), kickrank)];
            u = u*(rv.');
        end;
        r = size(v,2);        
        u = reshape(u, rx(i), n(i), r);
        v = reshape(v.', r, n(i+1), rx(i+2));
        
        % Recompute phi. Here are right phis
        phia{i+1} = compute_next_Phi(phia{i+2}, v, crA{i+1}, v, 'rl');
        phiy{i+1} = compute_next_Phi(phiy{i+2}, v, [], cry{i+1}, 'rl');        
    end;
    
    % Stuff back
    rx(i+1) = r;
    crx{i} = u;
    crx{i+1} = v;
    
    if (verb>2)
        somedata{1}(i,(swp-1)*2+1.5-dir/2) = toc(t_dmrg_solve);
        if (verb>3)||(i==d-1) % each microstep is returned only if really asked for
            % Otherwise the memory will blow up
            x = cell2core(tt_tensor, crx);
            somedata{2}{i, (swp-1)*2+1.5-dir/2}=x;
        end;
    end;
    
    i = i+dir;
    
    % Reversing, residue check, etc
    cur_order(order_index) = cur_order(order_index) - dir;
    % New direction
    if (cur_order(order_index)==0)
        order_index = order_index+1;
               
        if (verb>0)
            fprintf('=dmrg_solve3= sweep %d{%d}, max_dx: %3.3e, max_res: %3.3e, max_iter: %d, mrank: %d\n', swp, order_index-1, max_dx, max_res, max_iter, max(rx));
        end;        
        
        if (last_sweep)
            break;
        end;
        
        %residue
        if (strcmp(trunc_norm, 'fro'))
            if (max_dx<tol)&&(verb<3)
                last_sweep=true;
            end;
        else
            if (max_res<tol)&&(verb<3)
                last_sweep=true; % comment out to test
            end;
        end;
        
        % Reset the error indicators after the sweep
        if ((dir+dirfilter)==0)||(dirfilter==0)
            max_res = 0;
            max_dx = 0;
            max_iter = 0;
            dx_old = dx;
        end;
        
        if (order_index>numel(cur_order)) % New global sweep
            cur_order = block_order;
            order_index = 1;
            if (last_sweep)
                cur_order = d-1;
            end;
            swp = swp+1;            
        end;
        
        dir = sign(cur_order(order_index));
        i = i+dir;
    end;
end;


x = cell2core(tt_tensor, crx);

end


function [Phi] = compute_next_Phi(Phi_prev, x, A, y, direction)
% Performs the recurrent Phi (or Psi) matrix computation
% Phi = Phi_prev * (x'Ay).
% If direction is 'lr', computes Psi
% if direction is 'rl', computes Phi
% A can be empty, then only x'y is computed.

if (strcmp(direction, 'rl'))
  % Revert ranks to perform the right-to-left recursion
  x = permute(x, [3, 2, 1]);
  y = permute(y, [3, 2, 1]);
  if (~isempty(A))
    A = permute(A, [4, 2, 3, 1]);
  end
end

rx1 = size(x,1); n = size(x,2); rx2 = size(x,3);
ry1 = size(y,1); m = size(y,2); ry2 = size(y,3);
if (~isempty(A))
  ra1 = size(A,1); ra2 = size(A,4);
else
  ra1 = 1; ra2 = 1;
end

Phi = reshape(Phi_prev, [rx1*ra1, ry1]);
y = reshape(y, [ry1, m*ry2]);
Phi = Phi*y;	% complexity §\mcommentfont$\mathcal{O}(n  r_x r_A r_y^2)$§
Phi = reshape(Phi, [rx1, ra1, m, ry2]);
Phi = permute(Phi, [2, 3, 1, 4]);
if (~isempty(A))
  Phi = reshape(Phi, [ra1*m, rx1*ry2]);
  A = permute(A, [4, 2, 1, 3]);
  A = reshape(A, [ra2*n, ra1*m]);
  Phi = A*Phi;	% complexity §\mcommentfont$\mathcal{O}(n^2  r_x r_A^2 r_y)$§
  Phi = reshape(Phi, [ra2, n, rx1, ry2]);
end
Phi = permute(Phi, [3, 2, 1, 4]);
Phi = reshape(Phi, [rx1*n, ra2*ry2]);
x = reshape(x, [rx1*n, rx2]);
Phi = (x')*Phi;	% complexity §\mcommentfont$\mathcal{O}(n  r_x^2 r_A r_y)$§
if (~isempty(A))
  Phi = reshape(Phi, [rx2, ra2, ry2]);
end
end


function [y]=bfun3(Phi1,B1,B2,Phi2, x)
% Computes (Phi1 * B1 * B2 * Phi2)*x
% Phi1 is of sizes ry1, rB1, rx1
% B1 is of sizes rB1, k1, m1, rB2
% B2 is of sizes rB2, k2, m2, rB3
% Phi2 is of sizes ry3, rB3, rx3
ry1 = size(Phi1,1); ry3 = size(Phi2,1);
rx1 = size(Phi1,3); rx3 = size(Phi2,3);
rb1=size(B1,1); rb2=size(B1,4); rb3 = size(B2, 4);
m1 = size(B1,3); m2 = size(B2,3);
k1 = size(B1,2); k2 = size(B2,2);

y = reshape(x, rx1, m1*m2*rx3);
Phi1 = reshape(Phi1, ry1*rb1, rx1);
y = Phi1*y; % size ry1*rb1,m1*m2*rx3 % cplx rb*rx^3*m^2
y = reshape(y, ry1, rb1*m1, m2, rx3);
y = permute(y, [2, 1, 3, 4]);
y = reshape(y, rb1*m1, ry1*m2*rx3);
B1 = permute(B1, [2, 4, 1, 3]);
B1 = reshape(B1, k1*rb2, rb1*m1);
y = B1*y; % size k1*rb2, ry1*m2*rx3 % cplx rb^2*rx^2*n^3
y = reshape(y, k1, rb2, ry1, m2, rx3);
y = permute(y, [2, 4, 3, 1, 5]);
y = reshape(y, rb2*m2, ry1*k1*rx3);
B2 = permute(B2, [2, 4, 1, 3]);
B2 = reshape(B2, k2*rb3, rb2*m2);
y = B2*y; % size k2*rb3, ry1*k1*rx3 % cplx rb^2*rx^2*n^3
y = reshape(y, k2, rb3, ry1*k1, rx3);
y = permute(y, [2, 4, 3, 1]);
y = reshape(y, rb3*rx3, ry1*k1*k2);
Phi2 = reshape(Phi2, ry3, rb3*rx3);
y = Phi2*y; % size ry3, ry1*k1*k2 % cplx rb*rx^3*n^2
y = y.';
y = reshape(y, ry1*k1*k2*ry3, 1);
end
