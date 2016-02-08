function [x,iter,resids] = fgmres(A, b, tol, varargin)
% Flexible GMRES method
%   [x,iter] = fgmres(A, b, tol, varargin)
%
% Tries to solve a linear system Ax=b via the flexible GMRES
%
% A should be given as a function handle @(x,tol)A(x), which computes the
% Matrix-Vector product A*x with accuracy tol.
%
% A and b can be passed as regular double arrays
%
% varargin may contain a sequence of tuning parameters of the form 
% 'parameter1_name', parameter1_value, 'parameter2_name', parameter2_value
% and so on. Available parameters and [default values] are the following.
%   'relaxation' [0]: Inexact Krylov relaxation power. The accuracy for 
%                       Krylov vectors is set as tol/(|r|/|b|)^relaxation, 
%                       where r is the current residual
%   'max_iters' [2]: Number of outer iterations
%   'restart' [20]: Number of inner iterations
%   'x0' [all-zeros vector]: Initial guess
%   'verb' [2]:  Verbosity. 0 - Silent, 1 - outer iteration info, 2 - inner
%                iteration info, 3 - iterations + return all solutions in td{2}
%   'tol_exit' [tol]:  Stopping tolerance: exit if |r|/|b|<tol_exit
%   'P' []: Preconditioning procedure in the same form as A, i.e. @(x,tol)P(x)
%            
%
% Returns the solution x iteration numbers iter and residuals resids
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
relaxation = 0;
% Number of outer iterations
max_iters = 2;
% Number of inner iterations
restart = 20;
% Initial guess
x = [];
% Verbosity
verb = 2;
% Stopping tolerance
tol_exit = tol;
% Prec 
P = [];

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
        case 'tol_exit'
            tol_exit = varargin{i+1};            
        case 'p'
            P = varargin{i+1};
        otherwise
            error('Unknown tuning parameter "%s"', varargin{i});
    end;
end;


% If we are given a single matrix, just use exact MV
if (isa(A, 'double'))
    A = @(x,tol)(A*x);
end;

if (~isempty(P))&&(isa(P, 'double'))
    P = @(x,tol)(P*x);
end;

n = size(b,1);

if (restart>n)
    restart = n;
end;

% Norm of the RHS
beta0 = norm(b);
if (beta0==0)
    x = zeros(n,1);
    iter = 0;
    resids = 0;
    return;
end;

% Check for the initial guess
if (isempty(x))
    x = zeros(n,1);
end;

t_gmres_start = tic;

% % Hessinberg matrix
% H = zeros(restart+1, restart);
% Householder data
W = zeros(restart+1,1);
R = zeros(restart, restart);
J = zeros(2,restart);

% Krylov subspace
V = cell(restart, 1);
% Preconditioned subspace
if (~isempty(P))
    Z = cell(restart, 1);
end;

if (nargout>2)
    resids = zeros(restart, max_iters);
end;

for it=1:max_iters
    % Krylov tolerance
    tol_kryl = tol;
    % Compute the initial residual
    if (it>1)||(norm(x)>0)
        Ax = A(x, tol);
        u = b - Ax;
    else
        u = b;
    end;
    beta = norm(u);
    if (verb>0)
        fprintf('==fgmres== iter=%d, |r| = %g, time=%g  \n', it-1, beta, toc(t_gmres_start));
    end;
    % Prepare the first Householder vector
    if (u(1)<0)
        beta = -beta;
    end;
    u(1) = u(1)+beta;    
    u = u/norm(u);
    % The first Hh entry
    W(1) = -beta;
    V{1} = u;
    
    for j=1:restart
        % Construct the last vector from the HH storage
        %  Form P1*P2*P3...Pj*ej.
        %  v = Pj*ej = ej - 2*u*u'*ej
        v = -2*conj(u(j))*u;
        v(j) = v(j) + 1;
        %  v = P1*P2*...Pjm1*(Pj*ej)
        for i = (j-1):-1:1
            v = v - 2*V{i}*(V{i}'*v);
        end
        %  Explicitly normalize v to reduce the effects of round-off.
        v = v/norm(v);
        
        if (isempty(P))
            w = A(v, tol_kryl);
        else
            % Keep the PrecVec separately
            Z{j} = P(v, tol_kryl);
            w = A(Z{j}, tol_kryl);
        end;
        
        % Orthogonalize the Krylov vector
        %  Form Pj*Pj-1*...P1*Av.
        for i = 1:j
            w = w - 2*V{i}*(V{i}'*w);
        end        
        
        % Update the rotators
        %  Determine Pj+1.
        if (j ~= length(w))
            %  Construct u for Householder reflector Pj+1.
            u = [zeros(j,1); w(j+1:end)];
            alpha = norm(u);
            if (alpha ~= 0)
                if (w(j+1)<0)
                    alpha = -alpha;
                end;
                u(j+1) = u(j+1) + alpha;
                u = u / norm(u);
                V{j+1} = u;
                
                %  Apply Pj+1 to v.
                %  v = v - 2*u*(u'*v);
                w(j+2:end) = 0;
                w(j+1) = -alpha;
            end
        end
        
        %  Apply Given's rotations to the newly formed v.
        for colJ = 1:j-1
            tmpv = w(colJ);
            w(colJ)   = conj(J(1,colJ))*w(colJ) + conj(J(2,colJ))*w(colJ+1);
            w(colJ+1) = -J(2,colJ)*tmpv + J(1,colJ)*w(colJ+1);
        end
        
        %  Compute Given's rotation Jm.
        if (j~=length(w))
            rho = norm(w(j:j+1));
            J(:,j) = w(j:j+1)./rho;
            W(j+1) = -J(2,j).*W(j);
            W(j) = conj(J(1,j)).*W(j);
            w(j) = rho;
            w(j+1) = 0;
        end
        
        R(:,j) = w(1:restart);    
        
        % Local and global residuals
        err = abs(W(j+1))/abs(beta);
        resid = abs(W(j+1))/beta0;
        % Relax the Krylov tolerance
        tol_kryl = tol/(err^relaxation);        
                    
        % report the residual
        if (nargout>2)
            resids(j,it)=resid;
        end;
                    
        if (verb>1)
            fprintf('==fgmres== iter=[%d,%d], resid=%3.3e, locerr=%3.3e, time=%g\n', it, j, resid, err, toc(t_gmres_start));
        end;

        if (resid<tol_exit); break; end;
    end;
    
    % Correction
    y = R(1:j,1:j) \ W(1:j);
    if (isempty(P))
        dx = -2*y(j)*V{j}*conj(V{j}(j));
        dx(j) = dx(j) + y(j);
        for i = j-1:-1:1
            dx(i) = dx(i) + y(i);
            dx = dx - 2*V{i}*(V{i}'*dx);
        end
    else
        dx = zeros(n,1);
        for i=j:-1:1
            dx = dx + y(i)*Z{i};
        end;
    end;
    x = x+dx;
    
    if (resid<tol_exit); break; end;
end;

if (verb>0)
    fprintf('==fgmres== iter=%d, resid=%3.3e, time=%g\n', it, resid, toc(t_gmres_start));
end;

% Output
if (nargout>1)
    iter = [it, j];
end;
if (nargout>2)
    resids = resids(:, 1:it);
    resids(j+1:restart,it)=resid;
end;

end

