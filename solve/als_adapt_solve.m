function [x]=als_adapt_solve(A, y, tol, varargin)
%Solution of linear systems in TT-format via DMRG iteration
%   [X,SWEEPS]=ALS_ADAPT_SOLVE(A,Y,TOL,OPTIONS) Attempts to solve the linear
%   system A*X = Y with accuracy/residual TOL using the two-sided ALS+KICK iteration.
%   Matrix A has to be given in the TT-format, right-hand side Y should be
%   given in the TT-format also. Options are provided in form
%   'PropertyName1',PropertyValue1,'PropertyName2',PropertyValue2 and so
%   on. The parameters are set to default (in brackets in the following) 
%   The list of option names and default values are:
%       o x0 - initial approximation [random rank-2 tensor] 
%       o P - preconditioner  [I]
%       o nswp - maximal number of DMRG sweeps [10]
%       o rmax - maximal TT-rank of the solution [1000]
%       o verb - verbosity level, 0-silent, 1-sweep info, 2-block info [1]
%       o max_full_size - maximal size of the local matrix to full solver 
%       [2500]
%       o local_prec: Local preconditioner, 'als' - ALS-Richardson
%       iteration, 'selfprec' (Saad selfpreconditioner) ['als']
%       o prec_compr - compression for local precs [1e-3]
%       o prec_tol - tolerance for local precs [1e-1]
%       o prec_iters - number of local iterations [15]
%       o use_self_prec - Use self_prec [ true | {false} ]
%       o gmres_iters - number of local gmres restarts [2]
%       o nrestart - dimension of local gmres [40]
%       Example:
%           d=8; f=8; 
%           mat=tt_qlaplace_dd(d*ones(1,f)); %Laplace in the QTT-format
%           rhs=tt_ones(2,d*f); Right-hand side of all ones
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

als_tol_high = 10;
als_tol_low = 4;
als_iters = 3;

resid_damp = 2; % Truncation error to true residual treshold
bot_conv = 0.1; % bottom convergence factor - if better, we can decrease dpow and drank
top_conv = 0.99; % top convergence factor - if worse, we have to increase dpow and drank


nswp=10;
local_restart=40;
local_iters=2;

local_prec = '';
local_prec_char = 0;
% local_prec = 'jacobi';

rmax=1000;
trunc_norm = 'residual';
trunc_norm_char = 1;
% trunc_norm = 'fro';

local_solver = 'gmres';
ismex = true;
% local_solver = 'pcg';

verb=1;
kickrank = 5;
x=[];
block_order = [];
Au = [];

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
        case 'local_prec'
            local_prec=varargin{i+1};
        case 'local_restart'
            local_restart=varargin{i+1};
        case 'local_iters'
            local_iters=varargin{i+1};
        case 'local_solver'
            local_solver=varargin{i+1};            
        case 'kickrank'
            kickrank=varargin{i+1};
        case  'max_full_size'
            max_full_size=varargin{i+1};
        case 'step_dpow'
            step_dpow=varargin{i+1};
        case 'min_dpow'
            min_dpow=varargin{i+1};
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
        case 'als_tol_low'
            als_tol_low=varargin{i+1};
        case 'als_iters'
            als_iters=varargin{i+1};              
            
        otherwise
            error('Unrecognized option: %s\n',varargin{i});
    end
end

if (strcmp(local_prec, 'cjacobi')); local_prec_char = 1;  end;
if (strcmp(trunc_norm, 'fro')); trunc_norm_char = 0; end;

tol2 = tol;

if (A.n~=A.m)
    error(' DMRG does not know how to solve rectangular systems!\n Use dmrg_solve3(ctranspose(A)*A, ctranspose(A)*f, tol) instead.');
end;

d = y.d;
n = A.n;
if (isempty(x))
    x = tt_rand(n, A.d, kickrank);
end;

if (isempty(block_order))
    block_order = [+(d), -(d)];
end;

if (norm(y)==0) % Solution is a ground state. Keep it normalized
    x = x/norm(x);
end;

ry = y.r;
ra = A.r;
rx = x.r;

cry = core2cell(y);
crA = core2cell(A);
crx = core2cell(x);

phia = cell(d+1,1); phia{1}=1; phia{d+1}=1;
phiy = cell(d+1,1); phiy{1}=1; phiy{d+1}=1;

% This is for checking the residual via averaging
% cphia = cell(d+1,1); cphia{1}=1; cphia{d+1}=1;
% cphiy = cell(d+1,1); cphiy{1}=1; cphiy{d+1}=1;


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
    
    % For residual-check
%     cphia{i} = compute_next_Phi(cphia{i+1}, ones(1, n(i)), crA{i}, cr, 'rl');
%     cphiy{i} = compute_next_Phi(cphiy{i+1}, ones(1, n(i)), [], cry{i}, 'rl');
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
    Phi1 = phia{i}; Phi2 = phia{i+1};
    A1 = crA{i};
    % RHS
    rhs = phiy{i};
    rhs = rhs*reshape(cry{i}, ry(i), n(i)*ry(i+1));
    rhs = reshape(rhs, rx(i)*n(i), ry(i+1));
    rhs = rhs*(phiy{i+1}.');
    rhs = reshape(rhs, rx(i)*n(i)*rx(i+1),1);
    norm_rhs = norm(rhs);
    % sol_prev
    sol_prev = reshape(crx{i}, rx(i)*n(i)*rx(i+1), 1);
    
%     real_tol = (tol/(d^dpows(i)))/resid_damp;
%     if (last_sweep)
        real_tol = (tol2/sqrt(d))/resid_damp;
%     end;
    
    if (rx(i)*n(i)*rx(i+1)<max_full_size) % Full solution
        %      |     |    |
        % B = Phi1 - A1 - Phi2
        %      |     |    |
        B = reshape(permute(Phi1, [1, 3, 2]), rx(i)*rx(i), ra(i));
        B = B*reshape(A1, ra(i), n(i)*n(i)*ra(i+1));
        B = reshape(B, rx(i), rx(i), n(i), n(i), ra(i+1));
        B = permute(B, [1, 3, 2, 4, 5]);
        B = reshape(B, rx(i)*n(i)*rx(i)*n(i), ra(i+1));
        B = B*reshape(permute(Phi2, [2, 1, 3]), ra(i+1), rx(i+1)*rx(i+1));
        B = reshape(B, rx(i)*n(i), rx(i)*n(i), rx(i+1), rx(i+1));
        B = permute(B, [1, 3, 2, 4]);
        B = reshape(B, rx(i)*n(i)*rx(i+1), rx(i)*n(i)*rx(i+1));
        if (norm_rhs~=0)
            res_prev = norm(B*sol_prev-rhs)/norm_rhs;
            if (res_prev>real_tol)
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
        else
            % rhs = 0: we are looking for a ground state
            res_prev = norm(B*sol_prev);
            if (res_prev>real_tol)
                B2 = B+eye(rx(i)*n(i)*rx(i+1))*mean(abs(Phi1(:)))*mean(abs(A1(:)))*mean(abs(Phi2(:)));
                sol_prev2 = sol_prev;
                for it=1:local_restart
                    sol = B2 \ sol_prev2;
                    sol = sol/norm(sol);
                    res_new = norm(B*sol);
                    if (strcmp(trunc_norm, 'fro'))
                        if (norm(sol-sol_prev2)<real_tol); break; end;
                    else
                        if (res_new<real_tol); break; end;    
                    end;
                    sol_prev2 = sol;
                end;
%                 [sol,L]=eig(B);
%                 [v,l]=min(abs(diag(L)));
%                 sol = sol(:,l);
                flg = 0;
                iter=it;
            else
                sol = sol_prev;
                res_new = res_prev;
                flg=0;
                iter=0;
            end;
        end;
        
    else % Structured solution.
        
        if (norm_rhs~=0)
        res_prev = norm(bfun3(Phi1, A1, Phi2, sol_prev) - rhs)/norm_rhs;
        
        if (res_prev>real_tol)
            if (~ismex)
            if (strcmp(local_prec, 'jacobi')||strcmp(local_prec, 'seidel')) %&&(mod(order_index, 2)==1)   % &&(~last_sweep)
                % Prepare the Jacobi prec - on maximal rank
                if (rx(i)>rx(i+1))
                    jacPhi1 = reshape(permute(Phi1, [1, 3, 2]), rx(i)*rx(i), ra(i));
                    ind = (1:rx(i)) + (0:rx(i)-1)*rx(i); % diagonal elements
                    jacPhi = jacPhi1(ind,:);
                    % The second block, unfortunately, to compute explicitly
                    jacB = reshape(A1, ra(i)*n(i)*n(i), ra(i+1));
                    jacB = jacB*reshape(permute(Phi2, [2, 1, 3]), ra(i+1), rx(i+1)*rx(i+1));
                    jacB = reshape(jacB, ra(i), n(i), n(i), rx(i+1), rx(i+1));
                    jacB = permute(jacB, [2, 4, 3, 5, 1]);
                    jacB = reshape(jacB, n(i)*rx(i+1), n(i)*rx(i+1), ra(i));
                    
                    % Lower and upper triangular parts for Gauss-Seidel
%                     gsU = Phi1;
%                     gsL = Phi1;
%                     for k=1:rx(i)
% %                         gsU(k:rx(i), :, k)=0;
%                         gsL(1:k, :, k)=0;
%                     end;
                    
                    
                    jacBlocks = cell(rx(i),1);
%                     jacBlocks=zeros(n(i)*n(i+1)*rx(i+2), n(i)*n(i+1)*rx(i+2), rx(i));
                    for k=1:rx(i)
                        jacBlocks{k}=zeros(n(i)*rx(i+1), n(i)*rx(i+1));
                    end;
                    % So stupid preallocation of sparsity
%                     jacP = 0*speye(rx(i)*n(i)*n(i+1)*rx(i+2));
                    % Sum over ra
                    for k=1:ra(i)
                        for m=1:rx(i)
                            jacBlocks{m} = jacBlocks{m} + jacB(:,:,k)*jacPhi(m,k);
                        end;
                    end;
                    
                    jacdir = -1;
%                     Stuff into the prec
                    for m=1:rx(i)
%                         em = [zeros(m-1,1); 1; zeros(rx(i)-m, 1)];
                        jacBlocks{m} = inv(jacBlocks{m});
%                         jacBlocks{m} = jacBlocks{m}.';
%                         jacP = jacP + kron(jacBlocks{m}, spdiags(em, 0, rx(i), rx(i)));
                    end;
                    
                else
                    
                    jacPhi2 = reshape(permute(Phi2, [2, 1, 3]), ra(i+1), rx(i+1)*rx(i+1));
                    ind = (1:rx(i+1)) + (0:rx(i+1)-1)*rx(i+1); % diagonal elements
                    jacPhi = jacPhi2(:,ind);
                    % The second block, unfortunately, to compute explicitly
                    jacB = reshape(permute(Phi1, [1, 3, 2]), rx(i)*rx(i), ra(i));
                    jacB = jacB*reshape(A1, ra(i), n(i)*n(i)*ra(i+1));
                    jacB = reshape(jacB, rx(i), rx(i), n(i), n(i), ra(i+1));
                    jacB = permute(jacB, [1, 3, 2, 4, 5]);
                    jacB = reshape(jacB, rx(i)*n(i), rx(i)*n(i), ra(i+1));
                    
                    % Lower and upper triangular parts for Gauss-Seidel
                    if (strcmp(local_prec, 'seidel'))
                        gsL = Phi2;
                        for k=1:rx(i+1)
                            gsL(1:k, :, k)=0;
                        end;
                    end;
                    
                    jacBlocks = cell(rx(i+1),1);
                    for k=1:rx(i+1)
                        jacBlocks{k}=zeros(rx(i)*n(i), rx(i)*n(i));
                    end;
                    % Sum over ra
                    for k=1:ra(i+1)
                        for m=1:rx(i+1)
                            jacBlocks{m} = jacBlocks{m} + jacB(:,:,k)*jacPhi(k,m);
                        end;
                    end;

                    
                    jacdir = 1;
                    % Stuff into the prec
                    for m=1:rx(i+1)
                        jacBlocks{m} = inv(jacBlocks{m});
                    end;
                end; % of rx(i)<rx(i+2)
                
            elseif (strcmp(local_prec, 'cjacobi')) %&&(mod(order_index, 2)==1)
                
                % Jacobi prec on phys mode only
                jacPhi1 = reshape(permute(Phi1, [1, 3, 2]), rx(i)*rx(i), ra(i));
                ind = (1:rx(i)) + (0:rx(i)-1)*rx(i); % diagonal elements
                jacPhi1 = jacPhi1(ind,:);
                jacPhi2 = reshape(permute(Phi2, [2, 1, 3]), ra(i+1), rx(i+1)*rx(i+1));
                ind = (1:rx(i+1)) + (0:rx(i+1)-1)*rx(i+1); % diagonal elements
                jacPhi2 = jacPhi2(:,ind);
                
%                 tic;
%                 jacBlocks = cell(rx(i), rx(i+1));
                jacBlocks = zeros(n(i),n(i),rx(i),rx(i+1));
%                 jacBlocks = sparse(rx(i)*n(i)*rx(i+1), rx(i)*n(i)*rx(i+1));
                for k2=1:rx(i+1)
                    for k1=1:rx(i)
                        curblock = jacPhi1(k1,:)*reshape(A1, ra(i), n(i)*n(i)*ra(i+1));
                        curblock = reshape(curblock, n(i)*n(i), ra(i+1));
                        curblock = curblock*jacPhi2(:,k2);
                        curblock = reshape(curblock, n(i), n(i));
                        curblock = inv(curblock);
                        jacBlocks(:,:,k1,k2)=curblock;
%                         ind = k1+((1:n(i))-1)*rx(i)+(k2-1)*rx(i)*n(i);
%                         jacBlocks(ind, ind) = curblock;
%                         jacBlocks{k1,k2} = jacPhi1(k1,:)*reshape(A1, ra(i), n(i)*n(i)*ra(i+1));
%                         jacBlocks{k1,k2} = reshape(jacBlocks{k1,k2}, n(i)*n(i), ra(i+1));
%                         jacBlocks{k1,k2} = jacBlocks{k1,k2}*jacPhi2(:,k2);
%                         jacBlocks{k1,k2} = reshape(jacBlocks{k1,k2}, n(i), n(i));
%                         jacBlocks{k1,k2} = inv(jacBlocks{k1,k2});
                    end;
                end;
%                 jacgentime = toc
            else
                jacBlocks = [];
            end;
            
            
            drhs = rhs - bfun3(Phi1, A1, Phi2, sol_prev);
            if (isempty(jacBlocks))
                mvfun = @(v)bfun3(Phi1, A1, Phi2, v);
            else
                if (strcmp(local_prec, 'seidel'))
                    mvfun = @(v)bfun3(Phi1, A1, Phi2, gsfun(jacBlocks, gsL, A1, Phi1, v,jacdir));
                elseif (strcmp(local_prec, 'jacobi'))
                    mvfun = @(v)bfun3(Phi1, A1, Phi2, jacfun(jacBlocks,v,jacdir));
                else
                    mvfun = @(v)bfun3(Phi1, A1, Phi2, cjacfun(jacBlocks,v,rx(i),rx(i+1)));
%                     mvfun = @(v)bfun3(Phi1, A1, Phi2, jacBlocks*v);
                end;
            end;
            
            % Run the iterative solution
%             tic;
            fprintf('HERE!');
            if (strcmp(local_solver, 'gmres'))
                [dsol,flg,RELRES,iter] = gmres(mvfun, drhs, local_restart, min(real_tol/res_prev,1), local_iters);


                iter = (iter(1)-1)*local_restart + iter(2);
            else
                [dsol,flg,RELRES,iter] = pcg(mvfun, drhs, min(real_tol/res_prev,1), local_iters*local_restart);
            end;
%             soltime = toc
            
            if (~isempty(jacBlocks))
                if (strcmp(local_prec, 'seidel'))
                    dsol = gsfun(jacBlocks, gsL, A1, Phi1,  dsol ,jacdir);
                elseif (strcmp(local_prec, 'jacobi'))
                    dsol = jacfun(jacBlocks,dsol,jacdir);
                else
                    dsol = cjacfun(jacBlocks,dsol,rx(i),rx(i+1));
%                     dsol = jacBlocks*dsol;
                end;
            end;
            
            sol = sol_prev + dsol;            
            
            else
                                
%             sol = solve3d(permute(Phi1,[1,3,2]), A1, permute(Phi2, [1,3,2]), rhs, real_tol, trunc_norm_char, sol_prev, local_prec_char, local_restart, local_iters, 1);
            sol = solve3d_2(permute(Phi1,[1,3,2]), A1, permute(Phi2, [3,2,1]), rhs, real_tol, trunc_norm_char, sol_prev, local_prec_char, local_restart, local_iters, 1);

            flg=0;
            iter = 0;
%             if (res_new>real_tol); flg=1; end;
            end;
            
            res_new = norm(bfun3(Phi1, A1, Phi2, sol) - rhs)/norm_rhs;
        else
            sol = sol_prev;
            res_new = res_prev;
            flg=0;
            iter=0;
        end;
        
        
        else % Ground state
            res_prev = norm(bfun3(Phi1, A1, Phi2, sol_prev));
            
            if (res_prev>real_tol)
                sol_prev2 = sol_prev;
                Phi1mex = zeros(rx(i), rx(i), ra(i)+1);
                Phi1mex(1:rx(i), 1:rx(i), 1:ra(i)) = permute(Phi1,[1,3,2]);
                Phi1mex(1:rx(i), 1:rx(i), ra(i)+1) = eye(rx(i))*mean(abs(Phi1(:)));
                Phi2mex = zeros(rx(i+1), rx(i+1), ra(i+1)+1);
                Phi2mex(1:rx(i+1), 1:rx(i+1), 1:ra(i+1)) = permute(Phi2,[1,3,2]);
                Phi2mex(1:rx(i+1), 1:rx(i+1), ra(i+1)+1) = eye(rx(i+1))*mean(abs(Phi2(:)));
                Phi2mex = permute(Phi2mex, [2,3,1]);
                A1mex = zeros(n(i),n(i), ra(i)+1, ra(i+1)+1);
                A1mex(1:n(i), 1:n(i), 1:ra(i), 1:ra(i+1)) = permute(A1, [2, 3, 1, 4]);
                A1mex(1:n(i), 1:n(i), ra(i)+1, ra(i+1)+1) = eye(n(i))*mean(abs(A1(:)));
                A1mex = permute(A1mex, [3, 1, 2, 4]);
                for it=1:local_iters
                    rhs = sol_prev2;
%                     sol = solve3d(Phi1mex, A1mex, Phi2mex, rhs, real_tol, trunc_norm_char, sol_prev2, local_prec_char, local_restart, local_iters, 0);
                    sol = solve3d_2(Phi1mex, A1mex, Phi2mex, rhs, real_tol, trunc_norm_char, sol_prev2, local_prec_char, local_restart, local_iters, 0);
                    sol = sol/norm(sol);
                    res_new = norm(bfun3(Phi1, A1, Phi2, sol));
                    if (strcmp(trunc_norm, 'fro'))
                        if (norm(sol-sol_prev2)<real_tol); break; end;
                    else
                        if (res_new<real_tol); break; end;
                    end;
                    sol_prev2 = sol;
                end;
                flg=0;               
                iter = it;
            else
                sol = sol_prev;
                res_new = res_prev;
                flg=0;
                iter=0;
            end;
        end;
    end;
    
    if (flg>0)
        fprintf('-warn- local solver did not converge at block %d\n', i);
    end;
    if (res_prev/res_new<resid_damp)&&(res_new>real_tol)&&(norm_rhs~=0)
        fprintf('--warn-- the residual damp was smaller than in the truncation\n');
    end;
    
    dx(i) = norm(sol-sol_prev)/norm(sol);
    max_dx = max(max_dx, dx(i));
    
    % The new core does not converge - increase rank
%     if (dx(i)/dx_old(i)>top_conv)&&(dx(i)>tol)
%         dranks(i)=dranks(i)+1;
%         dpows(i)=dpows(i)+step_dpow;
%     end;
%     % The new core converges well - try to decrease rank
%     if (dx(i)/dx_old(i)<bot_conv)||(dx(i)<tol)
%         dranks(i)=max(dranks(i)-1, 0);
%         dpows(i)=max(dpows(i)-step_dpow, min_dpow);
%     end;
%     
%     if (last_sweep)
%         dpows(i)=0.5;
%         dranks(i)=0;
%     end;
    
    % Check the residual
%     cPhi1 = cphia{i}; cPhi2 = cphia{i+1};
%     crhs = cphiy{i};
%     crhs = crhs*reshape(cry{i}, ry(i), n(i)*ry(i+1));
%     crhs = reshape(crhs, n(i), ry(i+1));
%     crhs = crhs*(cphiy{i+1}.');
%     cAsol = bfun3(cPhi1, A1, cPhi2, sol);
%     chkres = norm(cAsol-crhs)/norm(crhs);
    
    
%     chkres = res_prev;
    max_res = max(max_res, res_prev);
%     max_res = max(max_res, res_prev);
    max_iter = max(max_iter, iter);
    
    % Truncation
    if (dir>0) % left-to-right
        sol = reshape(sol, rx(i)*n(i), rx(i+1));
    else
        sol = reshape(sol, rx(i), n(i)*rx(i+1));
    end;
    
    if (norm_rhs==0)
        norm_rhs=1;
    end;
    
    if (kickrank>=0)
    [u,s,v]=svd(sol, 'econ');
    s = diag(s);
    
    if (strcmp(trunc_norm, 'fro')) % We are happy with L2 truncation (when? but let it be)
%         r = my_chop2(s, max(real_tol, res_new)*resid_damp*norm(s));        
        r = my_chop2(s, real_tol*resid_damp*norm(s));        
    else
        if (dir>0); r = min(rx(i)*n(i),rx(i+1));
        else r = min(rx(i), n(i)*rx(i+1)); end;        
        cursol = u(:,1:r)*diag(s(1:r))*(v(:,1:r)');
        if (rx(i)*n(i)*rx(i+1)<max_full_size)
            res = norm(B*cursol(:)-rhs)/norm_rhs;
        else
            res = norm(bfun3(Phi1, A1, Phi2, cursol)-rhs)/norm_rhs;
        end;
        bfuncnt = 1;
        if (res<max(real_tol, res_new)*resid_damp)
            drank = -1;
        else
            drank = 1;
        end;        
        while (r>0)&&(r<=numel(s))
            cursol = u(:,1:r)*diag(s(1:r))*(v(:,1:r)');
            if (rx(i)*n(i)*rx(i+1)<max_full_size)
                res = norm(B*cursol(:)-rhs)/norm_rhs;
            else
                res = norm(bfun3(Phi1, A1, Phi2, cursol)-rhs)/norm_rhs;
            end;
            bfuncnt=bfuncnt+1;
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
        
%         bfuncnt = 0;
%         % Residual trunc; First, bin-search
%         r1 = 1; r2 = numel(s); r = round((r1+r2)/2);
%         while (r2-r1>1)
%             cursol = u(:,1:r)*diag(s(1:r))*(v(:,1:r)');
%             if (rx(i)*n(i)*rx(i+1)<max_full_size)
%                 res = norm(B*cursol(:)-rhs)/norm_rhs;
%             else
%                 res = norm(bfun3(Phi1, A1, Phi2, cursol)-rhs)/norm_rhs;
%             end;
%             bfuncnt=bfuncnt+1;
%             if (res<max(real_tol, res_new)*resid_damp)
%                 r2 = r;
%             else
%                 r1 = r;
%             end;
%             r = round((r1+r2)/2);
%         end;
%         r = max(r-1,1);
%         % More accurate Linear search
%         while (r<=numel(s))
%             cursol = u(:,1:r)*diag(s(1:r))*(v(:,1:r)');
%             if (rx(i)*n(i)*rx(i+1)<max_full_size)
%                 res = norm(B*cursol(:)-rhs)/norm_rhs;
%             else
%                 res = norm(bfun3(Phi1, A1, Phi2, cursol)-rhs)/norm_rhs;
%             end;
%             bfuncnt=bfuncnt+1;
%             if (res<max(real_tol, res_new)*resid_damp)
%                 break;
%             end;
%             r = r+1;
%         end;
    end;
    
    % Artificial rank increasing
%     r = r+dranks(i);
    r = min(r, numel(s));
    r = min(r, rmax);
    
    else
        if (dir>0)
            [u,v]=qr(sol, 0);
            v=v';
            r = size(u,2);
            s = ones(r,1);
        else
            [v,u]=qr(sol.', 0);
            v=conj(v);
            u=u.';
            r = size(u,2);
            s = ones(r,1);            
        end;        
    end;
    
    if (verb>1)
        fprintf('=dmrg_solve3=   block %d{%d}, dx: %3.3e, res: %3.3e, bfuncnt: %d, r: %d\n', i, dir, dx(i), res_prev, bfuncnt, r);
    end;
        
    if (dir>0)&&(i<d) % left-to-right, kickrank, etc
        u = u(:,1:r);
        v = conj(v(:,1:r))*diag(s(1:r));
        
        if (kickrank>0)
        % Smarter kick: low-rank PCA in residual
        % Matrix: Phi1-A{i}, rhs: Phi1-y{i}, sizes rx(i)*n - ra(i+1)
        leftA = permute(Phi1, [1, 3, 2]);
        leftA = reshape(leftA, rx(i)*rx(i), ra(i));
        leftA = leftA*reshape(A1, ra(i), n(i)*n(i)*ra(i+1));
        leftA = reshape(leftA, rx(i), rx(i), n(i), n(i), ra(i+1));
        leftA = permute(leftA, [1, 3, 5, 2, 4]);
        leftA = reshape(leftA, rx(i)*n(i)*ra(i+1), rx(i)*n(i));
        lefty = phiy{i};
        lefty = lefty*reshape(cry{i}, ry(i), n(i)*ry(i+1));
        lefty = reshape(lefty, rx(i)*n(i), ry(i+1));
        leftresid = leftA*reshape(u*v', rx(i)*n(i), rx(i+1));
        leftresid = reshape(leftresid, rx(i)*n(i), ra(i+1)*rx(i+1));
        leftresid = [leftresid, -lefty];
        
%         leftresid = randn(rx(i)*n(i), kickrank);
        
        % The right rank is now ra*rx+rf; Extract kickrank PCAs
%         uk = zeros(rx(i)*n(i), min(kickrank, n(i)*rx(i)));
%         for j=1:min(kickrank, n(i)*rx(i))
%             uk(:,j) = uchol(leftresid.', 1);
%             uk(:,1:j) = reort(uk(:,1:(j-1)), uk(:,j));
%             leftresid = leftA*uk(:,j);
%             leftresid = reshape(leftresid, rx(i)*n(i), ra(i+1));
%         end;
        uk = uchol(leftresid.', kickrank+1);
        uk = uk(:,end:-1:max(end-kickrank+1,1));        
%         [uk, s2, v2]=svd(leftresid, 'econ');
%         u2'*uk
%         uk = uk(:, 1:min(kickrank, size(uk,2))); % <- this is our new basis, ha-ha
%         
%         % Residual ~ orthogonal complement to the current basis
%         proj = uk - u*((u')*uk);
%         res = norm(proj)/norm(uk);
%         max_res = max(max_res, res);
        
%         uk = randn(rx(i)*n(i), kickrank); % <- old random basis
        
        % kick
            [u,rv]=qr([u,uk], 0);
            radd = size(uk, 2);
            v = [v, zeros(rx(i+1), radd)];
            v = v*(rv.');
%             u = reort(u, uk);
        end;
        cr2 = crx{i+1};
        cr2 = reshape(cr2, rx(i+1), n(i+1)*rx(i+2));
        v = (v.')*cr2; % size r+radd, n2, r3
        
        r = size(u,2);
%         r = r+radd;
        
        u = reshape(u, rx(i), n(i), r);
        v = reshape(v, r, n(i+1), rx(i+2));
        
        % Recompute phi. Left ones, so permute them appropriately
        phia{i+1} = compute_next_Phi(phia{i}, u, crA{i}, u, 'lr');
        phiy{i+1} = compute_next_Phi(phiy{i}, u, [], cry{i}, 'lr');
        
        % residual-check
%         cphia{i+1} = compute_next_Phi(cphia{i}, ones(1,n(i)), crA{i}, u, 'lr');
%         cphiy{i+1} = compute_next_Phi(cphiy{i}, ones(1,n(i)), [], cry{i}, 'lr');
                
        % Stuff back
        rx(i+1) = r;
        crx{i} = u;
        crx{i+1} = v;
    elseif (dir<0)&&(i>1) % right-to-left
        u = u(:,1:r)*diag(s(1:r));
        v = conj(v(:,1:r));
        
        if (kickrank>0)
        % Smarter kick: low-rank PCA in residual
        % Matrix: Phi1-A{i}, rhs: Phi1-y{i}, sizes rx(i)*n - ra(i+1)
        rightA = permute(Phi2, [2, 1, 3]);
        rightA = reshape(rightA, ra(i+1), rx(i+1)*rx(i+1));
        rightA = reshape(A1, ra(i)*n(i)*n(i), ra(i+1))*rightA;
        rightA = reshape(rightA, ra(i), n(i), n(i), rx(i+1), rx(i+1));
        rightA = permute(rightA, [2, 4, 1, 3, 5]);
        rightA = reshape(rightA, n(i)*rx(i+1)*ra(i), n(i)*rx(i+1));
        righty = reshape(cry{i}, ry(i)*n(i), ry(i+1));
        righty = righty*(phiy{i+1}.');
        righty = reshape(righty, ry(i), n(i)*rx(i+1));
        rightresid = rightA*(reshape(u*v', rx(i), n(i)*rx(i+1)).');
        rightresid = reshape(rightresid, n(i)*rx(i+1), ra(i)*rx(i));
        rightresid = [rightresid, -(righty.')];
        
%         rightresid = randn(n(i)*rx(i+1), kickrank);
        
        % The right rank is now ra*rx+rf; Extract kickrank PCAs
%         uk = zeros(n(i)*rx(i+1), min(kickrank, n(i)*rx(i+1)));
%         for j=1:min(kickrank, n(i)*rx(i+1))
%             uk(:,j) = uchol(rightresid.', 1);
%             uk(:,1:j) = reort(uk(:,1:(j-1)), uk(:,j));
%             rightresid = rightA*uk(:,j);
%             rightresid = reshape(rightresid, n(i)*rx(i+1), ra(i));
%         end;
        
        
        uk = uchol(rightresid.', kickrank+1);
        uk = uk(:,end:-1:max(end-kickrank+1,1));
%         [uk, s2, v2]=svd(rightresid, 'econ');
%         u2'*uk
%         uk = uk(:, 1:min(kickrank, size(uk,2))); % <- this is our new basis, ha-ha        
        
%         % Residual ~ orthogonal complement to the current basis
%         proj = uk - v*((v')*uk);
%         res = norm(proj)/norm(uk);
%         max_res = max(max_res, res);
        
        
%         uk = randn(n(i)*rx(i+1), kickrank); % <- old random
        
        % kick
            [v,rv]=qr([v,uk], 0);
%             v = reort(v, uk);
            radd = size(uk, 2);
            u = [u, zeros(rx(i), radd)];
            u = u*(rv.');
        end;
        cr2 = crx{i-1};
        cr2 = reshape(cr2, rx(i-1)*n(i-1), rx(i));
        u = cr2*u;
        
        r = size(v,2);
%         r = r+radd;
        
        u = reshape(u, rx(i-1), n(i-1), r);
        v = reshape(v.', r, n(i), rx(i+1));
        
        % Recompute phi. Here are right phis
        phia{i} = compute_next_Phi(phia{i+1}, v, crA{i}, v, 'rl');
        phiy{i} = compute_next_Phi(phiy{i+1}, v, [], cry{i}, 'rl');        
        % Residual check
%         cphia{i} = compute_next_Phi(cphia{i+1}, ones(1,n(i)), crA{i}, v, 'rl');
%         cphiy{i} = compute_next_Phi(cphiy{i+1}, ones(1,n(i)), [], cry{i}, 'rl');        
        
        % Stuff back
        rx(i) = r;
        crx{i-1} = u;
        crx{i} = v;
    elseif ((dir>0)&&(i==d))||((dir<0)&&(i==1))
        % Just stuff back the last core
        sol = u(:,1:r)*diag(s(1:r))*(v(:,1:r)');
        sol = reshape(sol, rx(i), n(i), rx(i+1));
        crx{i} = sol;
    end;
    
    
    i = i+dir;
    
    % Reversing, residue check, etc
    cur_order(order_index) = cur_order(order_index) - dir;
    % New direction
    if (cur_order(order_index)==0)
        order_index = order_index+1;
        
%         max_res = norm(A*x-y)/norm(y);
        
%         Au = tt_mvk4(A, x, tol, 'y0', Au, 'verb', 0);
%         fprintf('=dmrg_solve3= sweep %d{%d}, real_res: %3.3e\n', swp, order_index-1, norm(Au-y)/norm(y));
        if (verb>0)
%             fprintf('=dmrg_solve3= sweep %d{%d}, max_dx: %3.3e, max_res: %3.3e, real_res: %3.3e, max_iter: %d, erank: %g\n', swp, order_index-1, max_dx, max_res, norm(Au-y)/norm(y), max_iter, sqrt(rx(1:d)'*(n.*rx(2:d+1))/sum(n)));
            fprintf('=dmrg_solve3= sweep %d{%d}, max_dx: %3.3e, max_res: %3.3e, max_iter: %d, erank: %g\n', swp, order_index-1, max_dx, max_res, max_iter, sqrt(rx(1:d)'*(n.*rx(2:d+1))/sum(n)));
        end;        
        
        if (last_sweep)
            break;
        end;
        
        if (kickrank<0)
            kickrank=kickrank-1;
        end;        

            if (strcmp(trunc_norm, 'fro'))
                if (max_dx<tol)&&(kickrank<=-als_iters)
                    last_sweep=true;
                    kickrank = 0;
                end;
                if (max_dx<tol*als_tol_low)&&(kickrank>0)
                    kickrank=-1;
                end;                
            else
                if (max_res<tol)&&(kickrank<=-als_iters)
                    kickrank = 0;
                    last_sweep=true;
                end;
%                 if (max_res<tol*als_tol_low)&&(kickrank>0)
                if (max_res<tol)&&(kickrank>0)
%                     kickrank=-1;
                    kickrank = 0;
                    last_sweep=true;
                end;
            end;        
        
            max_res = 0;
            max_dx = 0;
            max_iter = 0;               
%             dx_old = dx;
            
            
        if (order_index>numel(cur_order)) % New global sweep
            cur_order = block_order;
            order_index = 1;
            %residue
            if (last_sweep)
                cur_order = d-1;
            end;
            x = cell2core(x, crx);
            swp = swp+1;            
        end;
        
        dir = sign(cur_order(order_index));
        i = i+dir;
    end;
end;

x = cell2core(x, crx);

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


function [y]=bfun3(Phi1,B1,Phi2, x)
% Computes (Phi1 * B1 * Phi2)*x
% Phi1 is of sizes ry1, rB1, rx1
% B1 is of sizes rB1, k1, m1, rB2
% Phi2 is of sizes ry2, rB2, rx2
ry1 = size(Phi1,1); ry2 = size(Phi2,1);
rx1 = size(Phi1,3); rx2 = size(Phi2,3);
rb1=size(B1,1); rb2=size(B1,4); 
m1 = size(B1,3);
k1 = size(B1,2);

y = reshape(x, rx1, m1*rx2);
Phi1 = reshape(Phi1, ry1*rb1, rx1);
y = Phi1*y; % size ry1*rb1,m1*rx2 % cplx rb*rx^3*m^2
y = reshape(y, ry1, rb1*m1, rx2);
y = permute(y, [2, 1, 3]);
y = reshape(y, rb1*m1, ry1*rx2);
B1 = permute(B1, [2, 4, 1, 3]);
B1 = reshape(B1, k1*rb2, rb1*m1);
y = B1*y; % size k1*rb2, ry1*rx2 % cplx rb^2*rx^2*n^3
y = reshape(y, k1, rb2, ry1, rx2);
y = permute(y, [2, 4, 3, 1]);
y = reshape(y, rb2*rx2, ry1*k1);
Phi2 = reshape(Phi2, ry2, rb2*rx2);
y = Phi2*y; % size ry2, ry1*k1 % cplx rb*rx^3*n^2
y = y.';
y = reshape(y, ry1*k1*ry2, 1);
end


function [y] = jacfun(jacs, x, dir) % Jacobi prec
m = numel(jacs);
n = size(jacs{1},1);
if (dir==1)
    y = reshape(x, n, m);
    for i=1:m
        %     y(:,i) = jacs(:,:,i) \ y(:,i);
        y(:,i) = jacs{i} * y(:,i);
    end;
else
    y = reshape(x, m, n);
    y = y.';
    for i=1:m
        %     y(:,i) = jacs(:,:,i) \ y(:,i);
        y(:,i) = jacs{i}*y(:,i);
    end;
    y = y.';
end;
y = y(:);
end

function [y] = gsfun(jacs, gsL, A1, Phi_full, x, dir) % Gauss-seidel prec
m = numel(jacs);
n = size(jacs{1},1);
if (dir==1)
    % Left-to-right, split on r3
    y = reshape(x*0, n, m);
    rhs = x;
    rhs = reshape(rhs, n, m);
    
    for i=1:m
        yprev = bfun3(Phi_full, A1, gsL(i,:,1:i-1), y(:, 1:i-1)); % L*y^1
        y(:,i) = jacs{i}*(rhs(:,i) - yprev);
    end;   
else
    % Right-to-left, split on r1
    y = reshape(x*0, m, n);
    rhs = x;
    rhs = reshape(rhs, m, n);
        
    for i=1:m
        yprev = bfun3(gsL(i,:,1:i-1), A1, Phi_full, y(1:i-1, :)); % L*y^1
        y(i, :) = (rhs(i, :) - yprev)*jacs{i};
    end;      
end;
y = y(:);
end


function [y] = cjacfun(jacs, x,r1,r2)

% n = size(jacs{1},1);
n = size(jacs,1);
y = reshape(x, r1, n, r2);
y = permute(y, [2, 1, 3]);
% y = reshape(y, n, r1*r2);
% y = mat2cell(y, n, ones(1,r1*r2));
% for k=1:r1*r2
%     y{k} = jacs{k}*y{k};
% end;
% for k2=1:r2
%     for k1=1:r1
% %         y(:,k1,k2) = jacs{k1,k2}*y(:,k1,k2);
%         y(:,k1,k2) = jacs(:,:,k1,k2)*y(:,k1,k2);
%     end;
% end;

y = cjacmex(jacs, y(:));

% y = cell2mat(y);
y = reshape(y, n, r1, r2);
y = permute(y, [2, 1 ,3]);
y = y(:);
end
