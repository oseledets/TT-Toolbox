function [x,td] = tt_gmres(A, b, tol, maxout, maxin, eps_x, eps_z, M1, M2, M3, x0, verbose, varargin)
%TT-GMRES method
%   [x,td] = TT_GMRES(A, B, TOL, MAXOUT, MAXIN, EPS_X, 
%   EPS_Z, M1, M2, M3, X0, VERBOSE, VARARGIN)
%   GMRES in TT format. Solve a system A*x=b, up to the accuracy tol,
%   with maximum number of outer iterations maxout, the size of Krylov basis
%   maxin, compression of solution eps_x, compression of Krylov vectors
%   eps_z, optional LEFT preconditioner: [M1*[M2]*[M3]], initial guess [x0]
%   (default is 1e-308*b), [verbose] - if 1 or unspecified - print messages,
%   if 0 - silent mode, if 2 - print messages, and also log solutions to td.
%
%   Debug uutput td is currently syncronised with the amen and dmrg:
%       td{1} - cumulative CPU time,
%       td{2} - current solutions x,
%       td{3} - local residuals
%   If you want the older (also full gmres) syntax with local residuals as
%   the second output, please git checkout to some commit before 01.06.2014
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

[atype,afun,afcnstr] = tt_iterchk(A);

%%%%%%%%%%%%% parameters %%%%%%%%%%%%%%%%%%%
max_swp=12; % maximum sweeps in mvk2 allowed
max_sp_factor = 1.5; % maximum stagpoints (as a part of Krylov basis size) allowed
max_zrank_factor = 4; % maximum rank of Krylov vectors with respect to the rank of X.
derr_tol_for_sp = 1.0; % minimal jump of residual at one iteration, which is considered as
                       % a stagnation.
use_err_trunc = 1; % use relaxation, i.e. compression accuracy eps_z/err for Krylov vectors 
err_trunc_power = 1; % theoretical estimate - 1 - sometimes worse than something like 0.8
compute_real_res = 0; % Compute the real residual on each step (expensive!)
% extra param in tt_iterapp.m: matvec_alg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargout>1)
    td = cell(3,1);
    td{1} = zeros(maxin, maxout);
    td{2} = cell(maxin, maxout);
    td{3} = zeros(maxin, maxout);
end;

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

norm_f = tt_dot2(pre_b,pre_b)*0.5;
mod_norm_f = sign(norm_f)*mod(abs(norm_f), 10);
order_norm_f = norm_f - mod_norm_f;
cur_norm_f = exp(mod_norm_f);

if (existsx0)
    x = x0;
else
    x = core(tt_zeros(tt_size(b), max(size(b))));
end;

H = zeros(maxin+1, maxin);
v = cell(maxin, 1);

err=2;
max_err=Inf;
old_err = 1;
stagpoints=0;

for nitout=1:maxout
    max_rank=[];
    err_for_trunc=1;    
    Ax = tt_iterapp('mtimes',afun,atype,afcnstr,x,min(eps_z/err_for_trunc, 0.5), max_rank, max_swp,varargin{:});
    r = tt_axpy2(0,1, pre_b, 0,-1, Ax, min(eps_z/err_for_trunc, 0.5), max_rank);
    if (existsM1)
        r = tt_iterapp('mldivide',m1fun,m1type,m1fcnstr,r,min(eps_z/err_for_trunc, 0.5), max_rank, max_swp,varargin{:});
        if (existsM2)
            r = tt_iterapp('mldivide',m2fun,m2type,m2fcnstr,r,min(eps_z/err_for_trunc, 0.5), max_rank, max_swp,varargin{:});
        end;
        if (existsM3)
            r = tt_iterapp('mldivide',m3fun,m3type,m3fcnstr,r,min(eps_z/err_for_trunc, 0.5), max_rank, max_swp,varargin{:});
        end;
    end;        
    beta = tt_dot2(r,r)*0.5;
    mod_beta = sign(beta)*mod(abs(beta), 10);
    order_beta = beta - mod_beta;    
    cur_beta = exp(mod_beta);    
    if (verbose>0)
        real_beta = exp(beta);
        fprintf('   cur_beta = %g, real_beta = %g\n', cur_beta, real_beta);
    end;
%     if (nitout==1) 
%         cur_normb = cur_beta; 
%         order_normb = order_beta;
%     end;    
    v{1} = tt_scal2(r, -beta, 1);
    for j=1:maxin
        max_w_rank = 0;
        w = tt_iterapp('mtimes',afun,atype,afcnstr,v{j},min(eps_z/err_for_trunc, 0.5), max_rank, max_swp,varargin{:});                
        max_w_rank = max([max_w_rank; tt_ranks(w)]);
        if (existsM1)
            w = tt_iterapp('mldivide',m1fun,m1type,m1fcnstr,w,min(eps_z/err_for_trunc, 0.5), max_rank, max_swp,varargin{:});
            max_w_rank = max([max_w_rank; tt_ranks(w)]);
            if (existsM2)
                w = tt_iterapp('mldivide',m2fun,m2type,m2fcnstr,w,min(eps_z/err_for_trunc, 0.5), max_rank, max_swp,varargin{:});
                max_w_rank = max([max_w_rank; tt_ranks(w)]);
            end;
            if (existsM3)
                w = tt_iterapp('mldivide',m3fun,m3type,m3fcnstr,w,min(eps_z/err_for_trunc, 0.5), max_rank, max_swp,varargin{:});
                max_w_rank = max([max_w_rank; tt_ranks(w)]);
            end;            
        end;
        
        max_wrank = max(tt_ranks(w));
                
        for i=1:j
            H(i,j)=tt_dot(w, v{i});
%             w = tt_axpy2(0,1, w, log(abs(H(i,j))+1e-308), -1*sign(H(i,j)), v{i}, min(eps_z/err_for_trunc, 0.5), max_rank);
            max_w_rank = max([max_w_rank; tt_ranks(w)]);
            w = tt_axpy2(0,1, w, log(abs(H(i,j))+1e-308), -1*sign(H(i,j)), v{i}, eps_z, max_rank);
            % Orthogonality muss sein!
%             w = tt_axpy3(0,1, w, log(abs(H(i,j))+1e-308), -1*sign(H(i,j)), v{i});
        end;
       
        H(j+1,j) = sqrt(tt_dot(w, w));
        if (j<maxin)
            v{j+1}=tt_scal2(w, -log(H(j+1,j)), 1);
        end;        
        
        [UH,SH,VH]=svd(H(1:j+1, 1:j), 0);
        SH = diag(SH);
        sigma_min_H = min(SH);
        sigma_max_H = max(SH);
        if (verbose>0)       
            fprintf('   min(sigma(H)) = %g\n', sigma_min_H);
        end;
        SH(1:numel(find(SH>1e-100)))=1./SH(1:numel(find(SH>1e-100))); %pseudoinv(SH)
        SH = diag(SH);
        y = cur_beta*VH*SH*(UH(1,:)'); % y = pinv(H)*(beta*e1)

        err = log(norm(H(1:j+1, 1:j)*y-[cur_beta zeros(1,j)]', 'fro')/cur_norm_f+1e-308);
        err = err+order_beta-order_norm_f;        
        err = exp(err);
        if (use_err_trunc==1)
            err_for_trunc = log(norm(H(1:j+1, 1:j)*y-[cur_beta zeros(1,j)]', 'fro')/cur_beta+1e-308);
            err_for_trunc = exp(err_for_trunc);
            err_for_trunc = (err_for_trunc*maxin*sigma_max_H/sigma_min_H)^(err_trunc_power);
            err_for_trunc = min(err_for_trunc, 1);
        end;
        
        % report the residual
        if (nargout>1)
            td{3}(j,nitout)=err;
        end;
        
        x_new = x;
        max_x_rank = 0;
        for i=j:-1:1
            x_new = tt_axpy2(0,1, x_new, log(abs(y(i))+1e-308)+order_beta, sign(y(i)), v{i}, eps_x);
            max_x_rank = max([max_x_rank; tt_ranks(x_new)]);            
        end;
        % Report the current solution, if required.
        if (verbose>1)
            td{2}{j,nitout} = x_new;
        end;
        
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
        if (derr<derr_tol_for_sp); stagpoints=stagpoints+1; end;       
        if (verbose>0)
            if (compute_real_res==1)
                fprintf('iter = [%d,%d], derr = %3.2f, resid=%3.2e, real_res=%3.2e, rank_w=%d, rank_x=%d, sp=%d, time=%g\n', nitout, j, derr, err, res, max_wrank, max_xrank, stagpoints, toc(t0));
            else                
                fprintf('iter = [%d,%d], derr = %3.2f, resid=%3.2e, rank_w=%d, rank_x=%d, sp=%d, time=%g\n', nitout, j, derr, err, max_w_rank, max_x_rank, stagpoints, toc(t0));
            end;
        end;

        % Report time
    	td{1}(j,nitout) = toc(t0);

        if (err<max_err)
            x_good = x_new;
            max_err=err;
        end;
        if (err<tol); break; end;
        if (stagpoints>=maxin*max_sp_factor); break; end;
    end;
    x = x_good;

    if (err<tol); break; end;
    if (stagpoints>=maxin*max_sp_factor); break; end;
end;

if (nargout>1)
    td{1}=td{1}(:, 1:nitout);
    td{2}=td{2}(:, 1:nitout);
    td{3}=td{3}(:, 1:nitout);
    if (nitout==1)
        td{1}=td{1}(1:j);
        td{2}=td{2}(1:j);
        td{2}=td{3}(1:j);
    end;
end;

end
