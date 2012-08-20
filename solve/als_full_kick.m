function [x]=als_full_kick(A, y, tol, varargin)

% Inner parameters
max_full_size=2500;

% als_tol_high = 10;
als_tol_low = 1;
als_iters = 3;

resid_damp = 2; % Truncation error to true residual treshold

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


verb=1;
kickrank = 1;
x=[];
block_order = [];

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
        case 'kickrank'
            kickrank=varargin{i+1};
        case  'max_full_size'
            max_full_size=varargin{i+1};
        case 'resid_damp'
            resid_damp = varargin{i+1};
        case 'trunc_norm'
            trunc_norm = varargin{i+1};
        case 'block_order'
            block_order=varargin{i+1};
            
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
    
end;


last_sweep = false;
swp = 1;
i = 1;

dx = zeros(d,1);
max_res = 0;
max_dx = 0;
max_iter = 0;

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
    
    real_tol = (tol2/sqrt(d))/resid_damp;
    
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
                sol = solve3d(permute(Phi1,[1,3,2]), A1, permute(Phi2, [1,3,2]), rhs, real_tol, trunc_norm_char, sol_prev, local_prec_char, local_restart, local_iters, 1);
                
                flg=0;
                iter = 0;
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
                A1mex = zeros(n(i),n(i), ra(i)+1, ra(i+1)+1);
                A1mex(1:n(i), 1:n(i), 1:ra(i), 1:ra(i+1)) = permute(A1, [2, 3, 1, 4]);
                A1mex(1:n(i), 1:n(i), ra(i)+1, ra(i+1)+1) = eye(n(i))*mean(abs(A1(:)));
                A1mex = permute(A1mex, [3, 1, 2, 4]);
                for it=1:local_iters
                    rhs = sol_prev2;
                    sol = solve3d(Phi1mex, A1mex, Phi2mex, rhs, real_tol, trunc_norm_char, sol_prev2, local_prec_char, local_restart, local_iters, 0);
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
    
    max_res = max(max_res, res_prev);
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

        if (dir>0)
            r_prev = rx(i+1);
        else
            r_prev = rx(i);
        end;
        if strcmp(trunc_norm, 'fro'); res_new = 0; end;
        if (rx(i)*n(i)*rx(i+1)<max_full_size)
            [r, bfuncnt]=resid_chop(trunc_norm, max(real_tol, res_new)*resid_damp, u, s, v, rhs, B, [], [], r_prev);
        else
            [r, bfuncnt]=resid_chop(trunc_norm, max(real_tol, res_new)*resid_damp, u, s, v, rhs, A1, Phi1, Phi2, r_prev);
        end;

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
        x2 = crx{i+1};
        x2 = reshape(x2, rx(i+1), n(i+1)*rx(i+2));
        v = v.'*x2;
        v = reshape(v, r, n(i+1), rx(i+2));
        u = reshape(u, rx(i), n(i), r);
        % Stuff back
        rx(i+1) = r;
        crx{i} = u;
        crx{i+1} = v;
        
        if (kickrank>0)&&(mod(i,kickrank)==0)
            rold = size(crx{i+1});
            if (numel(rold)<3); rold(3)=1; end;
            x = cell2core(x, crx);
%             xfullold = x;
            resid = A*x-y; 	    
            resid = core2cell(resid);
            for j=d:-1:2
                cr = resid{j};
                r2 = size(cr,1); r3 = size(cr,3);
                cr = reshape(cr, r2, n(j)*r3);
                [cr,rv]=qr(cr.',0);
                cr2 = resid{j-1};
                r1 = size(cr2, 1);
                cr2 = reshape(cr2, r1*n(j-1), r2);
                cr2 = cr2*(rv.');
                r2 = size(cr, 2);
                resid{j} = reshape(cr.', r2, n(j), r3);
                resid{j-1} = reshape(cr2, r1, n(j-1), r2);
            end;
            for j=1:d-1
                cr = resid{j};
                r1 = size(cr,1); r2 = size(cr,3);
                cr = reshape(cr, r1*n(j), r2);
                [uk,sk,vk]=svd(cr,'econ');
                uk=uk(:,1);
                vk=sk(1,1)*(vk(:,1)');
%                uk = uchol(cr.', 1);
%                uk = uk(:,end);
%                cr = uk'*cr;
                cr2 = resid{j+1};
                r3 = size(cr2, 3);
                cr2 = reshape(cr2, r2, n(j+1)*r3);
                cr2 = vk*cr2;
                resid{j} = reshape(uk, r1, n(j), 1);
                resid{j+1} = reshape(cr2, 1, n(j+1), r3);
            end;
            resid = cell2core(tt_tensor, resid);
%             fprintf('resid: %3.3e\n', norm(resid)/norm(y));

%             resid1 = als_rank_one(x, 'a', A, 'f', y, 'maxit', 100);
%             resid2 = als_rank_one(x, 'a', A, 'f', y, 'maxit', 100);            
%             [norm(resid1)/norm(y), norm(A*x-y)/norm(y)]
%             norm(A*x-y-round(A*x-y, tol, 1))/norm(A*x-y)
%             x = x+round(resid1+resid2, tol, 1);
            x=x+resid;
            rx = x.r;
            crx = core2cell(x);
            if (i==d-1)
                crx{i+1}(rold(1)+1:rx(i+1), :) = 0;
%                 kickrank = 0;
            else
                crx{i+1}(rold(1)+1:rx(i+1), :, rold(3)+1:rx(i+2)) = 0;
            end;
%             keyboard;
            for j=1:i
                cr = crx{j};
                cr = reshape(cr, rx(j)*n(j), rx(j+1));
                [cr, rv]=qr(cr, 0);
                cr2 = crx{j+1};
                cr2 = reshape(cr2, rx(j+1), n(j+1)*rx(j+2));
                cr2 = rv*cr2;
                rx(j+1) = size(cr, 2);
                cr = reshape(cr, rx(j), n(j), rx(j+1));
                crx{j+1} = reshape(cr2, rx(j+1), n(j+1), rx(j+2));
                crx{j} = cr;
                phia{j+1} = compute_next_Phi(phia{j}, cr, crA{j}, cr, 'lr');
                phiy{j+1} = compute_next_Phi(phiy{j}, cr, [], cry{j}, 'lr');                
            end;
            for j=d:-1:i+2
                cr = crx{j};
                cr = reshape(cr, rx(j), n(j)*rx(j+1));
                [cr, rv]=qr(cr.', 0);
                cr2 = crx{j-1};
                cr2 = reshape(cr2, rx(j-1)*n(j-1), rx(j));
                cr2 = cr2*(rv.');
                rx(j) = size(cr, 2);
                cr = reshape(cr.', rx(j), n(j), rx(j+1));
                crx{j-1} = reshape(cr2, rx(j-1), n(j-1), rx(j));
                crx{j} = cr;
                phia{j} = compute_next_Phi(phia{j+1}, cr, crA{j}, cr, 'rl');
                phiy{j} = compute_next_Phi(phiy{j+1}, cr, [], cry{j}, 'rl');
            end;
%             x = cell2core(x, crx);
%             norm(x-xfullold)/norm(x)
%             keyboard;
        else
            % Recompute phi. Left ones, so permute them appropriately
            phia{i+1} = compute_next_Phi(phia{i}, u, crA{i}, u, 'lr');
            phiy{i+1} = compute_next_Phi(phiy{i}, u, [], cry{i}, 'lr');
        end;
    elseif (dir<0)&&(i>1) % right-to-left
        u = u(:,1:r)*diag(s(1:r));
        v = conj(v(:,1:r));
        
        x2 = crx{i-1};
        x2 = reshape(x2, rx(i-1)*n(i-1), rx(i));
        u = x2*u;
        v = reshape(v.', r, n(i), rx(i+1));
        u = reshape(u, rx(i-1), n(i-1), r);
        % Stuff back
        rx(i) = r;
        crx{i-1} = u;
        crx{i} = v;

%         if (kickrank>0)
%             x = cell2core(x, crx);
%             resid = als_rank_one(x, 'a', A, 'f', y);
%             x = x+resid;
%             rx = x.r;
%             crx = core2cell(x);
%             for j=1:i-2
%                 cr = crx{j};
%                 cr = reshape(cr, rx(j)*n(j), rx(j+1));
%                 [cr, rv]=qr(cr, 0);
%                 cr2 = crx{j+1};
%                 cr2 = reshape(cr2, rx(j+1), n(j+1)*rx(j+2));
%                 cr2 = rv*cr2;
%                 rx(j+1) = size(cr, 2);
%                 cr = reshape(cr, rx(j), n(j), rx(j+1));
%                 crx{j+1} = reshape(cr2, rx(j+1), n(j+1), rx(j+2));
%                 crx{j} = cr;
%                 phia{j+1} = compute_next_Phi(phia{j}, cr, crA{j}, cr, 'lr');
%                 phiy{j+1} = compute_next_Phi(phiy{j}, cr, [], cry{j}, 'lr');                
%             end;
%             for j=d:-1:i
%                 cr = crx{j};
%                 cr = reshape(cr, rx(j), n(j)*rx(j+1));
%                 [cr, rv]=qr(cr.', 0);
%                 cr2 = crx{j-1};
%                 cr2 = reshape(cr2, rx(j-1)*n(j-1), rx(j));
%                 cr2 = cr2*(rv.');
%                 rx(j) = size(cr, 2);
%                 cr = reshape(cr.', rx(j), n(j), rx(j+1));
%                 crx{j-1} = reshape(cr2, rx(j-1), n(j-1), rx(j));
%                 crx{j} = cr;
%                 phia{j} = compute_next_Phi(phia{j+1}, cr, crA{j}, cr, 'rl');
%                 phiy{j} = compute_next_Phi(phiy{j+1}, cr, [], cry{j}, 'rl');
%             end;          
%         else
            % Recompute phi. Here are right phis
            phia{i} = compute_next_Phi(phia{i+1}, v, crA{i}, v, 'rl');
            phiy{i} = compute_next_Phi(phiy{i+1}, v, [], cry{i}, 'rl');
%         end;
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
        
        if (verb>0)
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
	    if (max_res<tol)
		kickrank=0;
		last_sweep=true;
	    end;
%            if (max_res<tol)&&(kickrank<=-als_iters)
%                kickrank = 0;
%                last_sweep=true;
%            end;
%            if (max_res<tol*als_tol_low)&&(kickrank>0)
%                kickrank=-1;
%            end;
        end;
        
        max_res = 0;
        max_dx = 0;
        max_iter = 0;
        
        
        if (order_index>numel(cur_order)) % New global sweep
            cur_order = block_order;
            order_index = 1;
            %residue
            if (last_sweep)
                cur_order = d-1;
            end;
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


function [r, bfuncnt]=resid_chop(trunc_norm, tol, u, s, v, rhs, A1, Phi1, Phi2, r_prev)
if (strcmp(trunc_norm, 'fro')) % We are happy with L2 truncation (when? but let it be)
    r = my_chop2(s, tol*norm(s));
    bfuncnt = 0;
else
    norm_rhs = norm(rhs);
    r = min(r_prev, size(s,1));
    cursol = u(:,1:r)*diag(s(1:r))*(v(:,1:r)');
    if (isempty(Phi1))
        res = norm(A1*cursol(:)-rhs)/norm_rhs;
    else
        res = norm(bfun3(Phi1, A1, Phi2, cursol)-rhs)/norm_rhs;
    end;
    bfuncnt = 1;
    if (res<tol)
        drank = -1;
    else
        drank = 1;
    end;
    while (r>0)&&(r<=numel(s))
        cursol = u(:,1:r)*diag(s(1:r))*(v(:,1:r)');
        if (isempty(Phi1))
            res = norm(A1*cursol(:)-rhs)/norm_rhs;
        else
            res = norm(bfun3(Phi1, A1, Phi2, cursol)-rhs)/norm_rhs;
        end;
        bfuncnt=bfuncnt+1;
        if (drank>0)
            if (res<tol)
                break;
            end;
        else
            if (res>=tol)
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
end
