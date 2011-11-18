% function [x] = tt_gmres(A, b, tol, maxout, maxin, eps_x, eps_z, M1, M2, M3, x0, verbose, varargin)
% GMRES in TT format. Solve a system A*x=b, up to the accuracy tol,
% with maximum number of outer iterations maxout, the size of Krylov basis
% maxin, compression of solution eps_x, compression of Krylov vectors
% eps_z, optional LEFT preconditioner: [M1*[M2]*[M3]], initial guess [x0]
% (default is 1e-308*b), [verbose] - if 1 or unspecified - print messages,
% if 0 - silent mode

function [x,RESVEC,rw,rx] = tt_gmres(A, b, tol, maxout, maxin, eps_x, eps_z, M1, M2, M3, x0, verbose, varargin)
[atype,afun,afcnstr] = tt_iterchk(A);

%%%%%%%%%%%%% parameters %%%%%%%%%%%%%%%%%%%
max_swp=12; % maximum sweeps in mvk2 allowed
max_sp_factor = 0.5; % maximum stagpoints (as a part of Krylov basis size) allowed
max_zrank_factor = 4; % maximum rank of Krylov vectors with respect to the rank of X.
derr_tol_for_sp = 1.0; % minimal jump of residual at one iteration, which is considered as
                       % a stagnation.
use_err_trunc = 1; % use compression accuracy eps_z/err for Krylov vectors 
err_trunc_power = 0.8; % theoretical estimate - 1 - sometimes sucks
compute_real_res = 0; % Compute the real residual on each step (expensive!)
% extra param in tt_iterapp.m: matvec_alg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t0 = tic;

existsM1=1;
existsM2=1;
existsM3=1;
existsx0=1;
if ((nargin<8)||(isempty(M1)))
    existsM1=0;
end;
if ((nargin<9)||(isempty(M2)))
    existsM2=0;
end;
if ((nargin<10)||(isempty(M3)))
    existsM3=0;
end;
if ((nargin<11)||(isempty(x0)))
    existsx0=0;
end;

if ((nargin<12)||(isempty(verbose)))
    verbose=1;
end;

pre_b = b;
max_rank=[];

if (compute_real_res==1)
    pre_b2 = b;
end;

if (existsM1)
    [m1type,m1fun,m1fcnstr] = tt_iterchk(M1);
    if (compute_real_res==1)
    pre_b2 = tt_iterapp('mldivide',m1fun,m1type,m1fcnstr,pre_b2,min(eps_z, 0.5), max_rank, max_swp,varargin{:});
    end;
end;
if (existsM2)
    [m2type,m2fun,m2fcnstr] = tt_iterchk(M2);
    if (compute_real_res==1)
    pre_b2 = tt_iterapp('mldivide',m2fun,m2type,m2fcnstr,pre_b2,min(eps_z, 0.5), max_rank, max_swp,varargin{:});    
    end;
end;
if (existsM3)
    [m3type,m3fun,m3fcnstr] = tt_iterchk(M3);
    if (compute_real_res==1)
    pre_b2 = tt_iterapp('mldivide',m3fun,m3type,m3fcnstr,pre_b2,min(eps_z, 0.5), max_rank, max_swp,varargin{:});    
    end;
end;
% pre_b = tt_stabilize(pre_b,0);

% norm_f = tt_dot2(pre_b,pre_b)*0.5;
% mod_norm_f = sign(norm_f)*mod(abs(norm_f), 10);
% order_norm_f = norm_f - mod_norm_f;
% cur_norm_f = exp(mod_norm_f);

if (existsx0)
    x = x0;
else
%     x = pre_b;
%     x = tt_scal2(b, log(1e-308), 1);
    x = tt_zeros(max(size(b)), tt_size(b));
end;

% x_old = x;

H = zeros(maxin+1, maxin);
v = cell(maxin, 1);

if (nargout>1)
    RESVEC = ones(maxout, maxin);
end;
if (nargout>2)
    rw = zeros(maxout,maxin);
end;
if (nargout>2)
    rx = zeros(maxout,maxin);
end;
% z = cell(maxin, 1);

err=2;
max_err=2;
old_err = 1;
stagpoints=0;

for nitout=1:maxout
    max_rank=[];
    err_for_trunc=1;    
    Ax = tt_iterapp('mtimes',afun,atype,afcnstr,x,min(eps_z/err_for_trunc, 0.5), max_rank, max_swp,varargin{:});
    r = tt_axpy2(0,1, pre_b, 0,-1, Ax, min(eps_z/err_for_trunc, 0.5), max_rank);
    if (existsM1)
        r = tt_iterapp('mldivide',m1fun,m1type,m1fcnstr,r,min(eps_z/err_for_trunc, 0.5), max_rank, max_swp,varargin{:});
%         Ax = tt_stabilize(Ax,0);
        if (existsM2)
            r = tt_iterapp('mldivide',m2fun,m2type,m2fcnstr,r,min(eps_z/err_for_trunc, 0.5), max_rank, max_swp,varargin{:});
%             Ax = tt_stabilize(Ax,0);
        end;
        if (existsM3)
            r = tt_iterapp('mldivide',m3fun,m3type,m3fcnstr,r,min(eps_z/err_for_trunc, 0.5), max_rank, max_swp,varargin{:});
%             Ax = tt_stabilize(Ax,0);
        end;
    end;        
%     Ax = tt_stabilize(Ax,0);
%     r = tt_add(pre_b, tt_scal(Ax, -1));
%     r = tt_compr2(r, min(eps_z/err, 0.5), max_rank);
%     beta = sqrt(tt_dot(r, r))
    beta = tt_dot2(r,r)*0.5;
    mod_beta = sign(beta)*mod(abs(beta), 10);
    order_beta = beta - mod_beta;    
    cur_beta = exp(mod_beta);    
    if (verbose==1)
        real_beta = exp(beta);
        fprintf('   cur_beta = %g, real_beta = %g\n', cur_beta, real_beta);
    end;
%     if (beta<tol) break; end;
    if (nitout==1) 
        cur_normb = cur_beta; 
        order_normb = order_beta;
%         normb = beta;
    end;    
    v{1} = tt_scal2(r, -beta, 1);
%     v{1}=tt_stabilize(v{1},0);
    for j=1:maxin
        max_w_rank = 0;
        w = tt_iterapp('mtimes',afun,atype,afcnstr,v{j},min(eps_z/err_for_trunc, 0.5), max_rank, max_swp,varargin{:});                
        max_w_rank = max([max_w_rank; tt_ranks(w)]);
%         w=tt_stabilize(w,0);
%         max_Mx_rank = max(tt_ranks(w))
        if (existsM1)
            w = tt_iterapp('mldivide',m1fun,m1type,m1fcnstr,w,min(eps_z/err_for_trunc, 0.5), max_rank, max_swp,varargin{:});
            max_w_rank = max([max_w_rank; tt_ranks(w)]);
%             w=tt_stabilize(w,0);
            if (existsM2)
                w = tt_iterapp('mldivide',m2fun,m2type,m2fcnstr,w,min(eps_z/err_for_trunc, 0.5), max_rank, max_swp,varargin{:});
                max_w_rank = max([max_w_rank; tt_ranks(w)]);
%                 w=tt_stabilize(w,0);
%                 max_Px_rank = max(tt_ranks(w))                
            end;
            if (existsM3)
                w = tt_iterapp('mldivide',m3fun,m3type,m3fcnstr,w,min(eps_z/err_for_trunc, 0.5), max_rank, max_swp,varargin{:});
                max_w_rank = max([max_w_rank; tt_ranks(w)]);
%                 w=tt_stabilize(w,0);
            end;            
        end;
        
        max_wrank = max(tt_ranks(w));
                
%         wnew = w;
        for i=1:j
            H(i,j)=tt_dot(w, v{i});
%             w = tt_axpy2(0,1, w, log(abs(H(i,j))+1e-308), -1*sign(H(i,j)), v{i}, min(eps_z/err_for_trunc, 0.5), max_rank);
            max_w_rank = max([max_w_rank; tt_ranks(w)]);
            w = tt_axpy2(0,1, w, log(abs(H(i,j))+1e-308), -1*sign(H(i,j)), v{i}, eps_z, max_rank);
            % Orthogonality muss sein!
%             w = tt_axpy3(0,1, w, log(abs(H(i,j))+1e-308), -1*sign(H(i,j)), v{i});
        end;
%         w = tt_wround([], w, min(eps_z/err_for_trunc,0.5), 'rmax', max_rank, 'nswp', max_swp, 'verb', 1);
%         for i=1:j
%             wv = abs(tt_dot(w, v{i}));
%             if (wv>eps_z)
%                 fprintf('Warning: dot(w,v{%d}) = %3.3e\n', i, wv);
%             end;
%         end;
%         w=wnew;
        if (nargout>2)
            rw(nitout,j)=max_w_rank;
        end;
        
        H(j+1,j) = sqrt(tt_dot(w, w));
        if (j<maxin)
            v{j+1}=tt_scal2(w, -log(H(j+1,j)), 1);
%             v{j+1}=tt_stabilize(v{j+1},0);
        end;        
        
        [UH,SH,VH]=svd(H(1:j+1, 1:j), 0);
        SH = diag(SH);
        sigma_min_H = min(SH);
        sigma_max_H = max(SH);
        if (verbose==1)       
            fprintf('   min(sigma(H)) = %g\n', sigma_min_H);
        end;
        SH(1:numel(find(SH>1e-100)))=1./SH(1:numel(find(SH>1e-100))); %pseudoinv(SH)
        SH = diag(SH);
        y = cur_beta*VH*SH*(UH(1,:)'); % y = pinv(H)*(beta*e1)
%         y = beta*VH*SH*(UH(1,:)'); % y = pinv(H)*(beta*e1)

%         err = norm(H(1:j+1, 1:j)*y-[beta zeros(1,j)]', 'fro')/normb;
        err = log(norm(H(1:j+1, 1:j)*y-[cur_beta zeros(1,j)]', 'fro')/cur_normb+1e-308);
        err = err+order_beta-order_normb;
        err = exp(err);
        if (use_err_trunc==1)
%             err_for_trunc=err;
            err_for_trunc = log(norm(H(1:j+1, 1:j)*y-[cur_beta zeros(1,j)]', 'fro')/cur_beta+1e-308);
            err_for_trunc = exp(err_for_trunc);
            err_for_trunc = (err_for_trunc*maxin*sigma_max_H/sigma_min_H)^(err_trunc_power);
            err_for_trunc = min(err_for_trunc, 1);
%             err_for_trunc = err;
        end;
%         err_to_f = log(norm(H(1:j+1, 1:j)*y-[cur_beta zeros(1,j)]', 'fro')/cur_norm_f+1e-308);
%         err_to_f = err_to_f+order_beta-order_norm_f;
%         err_to_f = exp(err_to_f);        
        
        if (nargout>1)
            RESVEC(nitout,j)=err;
        end;
        
        
        x_new = x;
        max_x_rank = 0;
%         dx = tt_scal2(v{j}, log(abs(y(j))+1e-308)+order_beta, sign(y(j)));
        for i=1:j
%             dx = tt_add(dx, tt_scal2(v{i}, log(abs(y(i))+1e-308)+order_beta, sign(y(i))));
%             dx = tt_axpy2(0,1, dx, log(abs(y(i))+1e-308)+order_beta, sign(y(i)), v{i}, eps_x);
            x_new = tt_axpy2(0,1, x_new, log(abs(y(i))+1e-308)+order_beta, sign(y(i)), v{i}, eps_x);
%             x_new = tt_axpy3(0,1, x_new, log(abs(y(i))+1e-308)+order_beta, sign(y(i)), v{i});
            max_x_rank = max([max_x_rank; tt_ranks(x_new)]);            
        end;
        if (nargout>3)
            rx(nitout,j)=max_x_rank;
        end;
%         x_new = tt_wround([], x_new, eps_x, 'verb', 1);
%         x_new = tt_add(x, dx);
%         x_new = tt_compr2(x_new, eps_x);
%         x_new = tt_axpy2(0,1,x,0,1,dx,eps_x);
       
        max_xrank = max(tt_ranks(x_new));
        max_rank = max_zrank_factor*max_xrank;
        
        if (compute_real_res==1)
            res = tt_iterapp('mtimes',afun,atype,afcnstr,x_new,eps_z, max_rank, max_swp,varargin{:});
            if (existsM1)
                res = tt_iterapp('mldivide',m1fun,m1type,m1fcnstr,res,eps_z, max_rank, max_swp,varargin{:});
                if (existsM2)
                    res = tt_iterapp('mldivide',m2fun,m2type,m2fcnstr,res,eps_z, max_rank, max_swp,varargin{:});
                end;
                if (existsM3)
                    res = tt_iterapp('mldivide',m3fun,m3type,m3fcnstr,res,eps_z, max_rank, max_swp,varargin{:});
                end;
            end;
            res = tt_dist3(res, pre_b2)/sqrt(tt_dot(pre_b2,pre_b2));
        end;
        
        derr = old_err/err;
        old_err = err;
%         dx = tt_add(x_new, tt_scal(x_old, -1));
%         normx = sqrt(tt_dot(x_new,x_new));
%         x_old=x_new;
%         dx_norm = sqrt(tt_dot(dx, dx))/normx;
        if (derr<derr_tol_for_sp) stagpoints=stagpoints+1; end;       
        if (verbose==1)
%             err_for_trunc
            if (compute_real_res==1)
                fprintf('iter = [%d,%d], derr = %3.2f, resid=%3.2e, real_res=%3.2e, rank_w=%d, rank_x=%d, sp=%d, time=%g\n', nitout, j, derr, err, res, max_wrank, max_xrank, stagpoints, toc(t0));
            else                
                fprintf('iter = [%d,%d], derr = %3.2f, resid=%3.2e, rank_w=%d, rank_x=%d, sp=%d, time=%g\n', nitout, j, derr, err, max_w_rank, max_x_rank, stagpoints, toc(t0));
            end;
        end;
        if (err<max_err)
            x_good = x_new;
            max_err=err;
        end;
        if (err<tol) break; end;
        if (stagpoints>=maxin*max_sp_factor) break; end;
    end;
    x = x_good;

    if (err<tol) break; end;
    if (stagpoints>=maxin*max_sp_factor) break; end;
end;

if (nargout>1)
    RESVEC=RESVEC(1:nitout,:);
    if (nitout==1)
        RESVEC=RESVEC(:,1:j);
    end;
end;
if (nargout>2)
    rw=rw(1:nitout,:);
    if (nitout==1)
        rw=rw(:,1:j);
    end;
end;
if (nargout>3)
    rx=rx(1:nitout,:);
    if (nitout==1)
        rx=rx(:,1:j);
    end;
end;

end