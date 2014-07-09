function [x,td] = tt_gmres_block(A, b, tol, varargin)
%TT-GMRES method with a separated block storage of the solution.
%   [x,td] = tt_gmres_block(A, b, tol, varargin)
%
% Tries to solve a linear system Ax=b, where x and b are cell arrays of
% size 1 x M, each b{m} is a tt_tensor, such that the total vector would
% write as [b{1}; ... ; b{M}].
%
% A should be given as a function handle @(x,tol)A(x), which takes and return
% cell arrays of tt_tensors of the form given above, and computes the
% Matrix-Vector product with the approximation tolerance tol. Note that the 
% left preconditioner should be also incorporated into A, and hence in b.
%
% A and b can be passed as regular tt_matrix and tt_tensor, resp. However,
% the result will still be returned in the form cell(1,1).
%
% varargin may contain a sequence of tuning parameters of the form 
% 'parameter1_name', parameter1_value, 'parameter2_name', parameter2_value
% and so on. Available parameters and [default values] are the following.
%   'relaxation' [0.5]: Inexact Krylov relaxation power. The accuracy for 
%                       Krylov vectors is set as tol/(|r|/|b|)^relaxation, 
%                       where r is the current residual
%   'max_iters' [2]: Number of outer iterations
%   'restart' [20]: Number of inner iterations
%   'x0' [all-zeros vector]: Initial guess
%   'verb' [2]:  Verbosity. 0 - Silent, 1 - outer iteration info, 2 - inner
%                iteration info, 3 - iterations + return all solutions in td{2}
%   'axpy' [round(y+ax)]: A function to compute approximation to y = y+ax.
%                         Must be a function of the form @(a,x,y,tol)axpy,
%                         where a is a scalar, x and y are cells of
%                         tt_tensors as shown above, tol is the accuracy.
%   'tol_exit' [tol]:  Stopping tolerance: exit if |r|/|b|<tol_exit
%            
%
%   Debug output td is currently syncronised with the amen and dmrg:
%       td{1} - cumulative CPU time,
%       td{2} - current solutions x,
%       td{3} - local residuals
%
%
% Development: Sergey Dolgov ( sergey.v.dolgov@gmail.com ),
% Max Planck Institute for Mathematics in the Sciences, Leipzig.
%
%TT-Toolbox 2.2, 2009-2012
%
%This is TT Toolbox, written by Ivan Oseledets et al.
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------


% Inexact Krylov relaxation power. The accuracy for Krylov vectors is set
% as tol/(|r|/|b|)^relaxation, where r is the current residual
relaxation = 0.5;
% Number of outer iterations
max_iters = 2;
% Number of inner iterations
restart = 20;
% Initial guess
x = [];
% Verbosity
verb = 2;
% aX+Y function
axpy = @(a,x,y,tol)roundaxpy(a,x,y,tol);
% Stopping tolerance
tol_exit = tol;

for i=1:2:length(varargin)-1
    switch lower(varargin{i})
        case 'relaxation'
            relaxation = varargin{i+1};
        case 'max_iters'
            max_iters = varargin{i+1};
        case 'restart'
            restart = varargin{i+1};
        case 'x0'
            x = varargin{i+1};
        case 'verb'
            verb = varargin{i+1};
        case 'axpy'
            axpy = varargin{i+1};
        case 'tol_exit'
            tol_exit = varargin{i+1};            
        otherwise
            error('Unknown tuning parameter "%s"', varargin{i});
    end;
end;

% Block size
if (isa(b, 'tt_tensor'))
    M = 1;
    r = {b};
else
    M = numel(b);
    r = b;
end;

% If we are given a single tt_matrix, just use exact MV+round
if (isa(A, 'tt_matrix'))
    A = @(x,tol){round(A*x{1}, tol)};
end;

% Debug info
if (nargout>1)
    td = cell(3,1);
    td{1} = zeros(restart, max_iters);
    td{2} = cell(restart, max_iters);
    td{3} = zeros(restart, max_iters);
end;

% Check for the initial guess
if (isempty(x))
    x = cell(1,M);
    for m=1:M
        x{m} = tt_zeros(r{m}.n);
    end;
else
    if (isa(x,'tt_tensor'))
        x = {x};
    end;
end;

t_gmres_start = tic;

nrm_b = zeros(1,M);
for m=1:M
    nrm_b(m) = norm(r{m});
end;

% Hessinberg matrix
H = zeros(restart+1, restart);
% Krylov subspace
V = cell(restart, 1);

max_x_rank = 0;
max_w_rank = 0;

for it=1:max_iters
    % Krylov tolerance
    tol_kryl = tol;
    % Compute the initial residual
    if (it>1)||(cellnorm(x)>0)
        Ax = A(x, tol);
        V{1} = axpy(-1, Ax, r, tol);
    else
        V{1} = r;
    end;
    beta = cellnorm(V{1});
    if (it==1); beta0 = beta; end;
    if (verb>0)
        fprintf('==tt_gmres== iter=%d, |r| = %g, rank_w=%d, rank_x=%d, time=%g  \n', it-1, beta, max_w_rank, max_x_rank, toc(t_gmres_start));
    end;
    for m=1:M
        V{1}{m} = V{1}{m}/beta;
    end;
    for j=1:restart
        w = A(V{j}, tol_kryl);
             
        % Modified Gram-Schmidt
        for i=1:j
            H(i,j)=celldot(V{i}, w);
            w = axpy(-H(i,j), V{i}, w, tol_kryl);
        end;
       
        % Store the next Krylov vector
        H(j+1,j) = cellnorm(w);
        if (j<restart)
            V{j+1} = w;
            for m=1:M
                V{j+1}{m} = V{j+1}{m}/H(j+1,j);
                max_w_rank = max(max_w_rank, max(V{j+1}{m}.r));
            end;
        end;
        
        % Solve the LS-system by a brute force QR
        [Q,R]=qr(H(1:j+1, 1:j), 0);
        y = R\(Q(1,:)'*beta);
        err = norm(H(1:j+1, 1:j)*y - [beta;zeros(j,1)])/beta;
        resid = norm(H(1:j+1, 1:j)*y - [beta;zeros(j,1)])/beta0;
        % Relax the Krylov tolerance
        tol_kryl = tol/(err^relaxation);
        
        % report the residual
        if (nargout>1)
            td{3}(j,it)=err;
        end;
        
        % Report the current solution, if required
        if (verb>2)
            % Correction
            x_new = x;
            for i=j:-1:1
                x_new = axpy(y(i), V{i}, x_new, tol);
            end;
            td{2}{j,it} = x_new;
            for m=1:M
                max_x_rank = max(max_x_rank, max(x_new{m}.r));
            end;
        end;
        % Report time
    	td{1}(j,it) = toc(t_gmres_start);
              
        if (verb>1)
            fprintf('==tt_gmres== iter=[%d,%d], resid=%3.3e, locerr=%3.3e, rank_w=%d, rank_x=%d, time=%g\n', it, j, resid, err, max_w_rank, max_x_rank, td{1}(j,it));
        end;

        if (resid<tol_exit); break; end;
    end;
    
    % Correction
    x_new = x;
    for i=j:-1:1
        x_new = axpy(y(i), V{i}, x_new, tol);
    end;
    for m=1:M
        max_x_rank = max(max_x_rank, max(x_new{m}.r));
    end;
    x = x_new;
    
    if (resid<tol_exit); break; end;
end;

if (nargout>1)
    td{1}=td{1}(:, 1:it);
    td{2}=td{2}(:, 1:it);
    td{3}=td{3}(:, 1:it);
    if (it==1)
        td{1}=td{1}(1:j);
        td{2}=td{2}(1:j);
        td{2}=td{3}(1:j);
    end;
end;

end

% Just the rounding for axpy
function y = roundaxpy(a,x,y,tol)
M = numel(x);
for m=1:M
    y{m} = a*x{m}+y{m};
    y{m} = round(y{m}, tol);
end;
end

% Compute the second norm of the cell of tt_tensors
function nrm = cellnorm(x)
nrm = 0;
M = numel(x);
for m=1:M
    nrm = nrm+norm(x{m})^2;
end;
nrm = sqrt(nrm);
end

% Compute the dot product of cells of tt_tensors
function sc = celldot(x,y)
sc = 0;
M = numel(x);
for m=1:M
    sc = sc+dot(x{m},y{m});
end;
end