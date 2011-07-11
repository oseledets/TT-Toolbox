function [u]=newton_solve(B, M, A, f, AMh, tol, u0)

maxit=20;
n = size(B,1);
if (nargin<7)||(isempty(u0))
    u0 = 1e-5*ones(n,1);
end;

u = u0;

for i=1:maxit
    Bu = B*u;
    Mu = M*u;
    Au = A*u;
    
    resid = Bu + Mu.*Au - f;
    err_old = norm(resid);
    
    diag_Au = spdiags(Au, 0, n, n);
    diag_Mu = spdiags(Mu, 0, n, n);
    
    grad_resid = B + diag_Mu*A + diag_Au*M; % [i,\hat{i}] - matrix
    Grad = (grad_resid')*resid;
    
    Hess = AMh*resid; % size \tilde{i}*\hat{i},1
    Hess = reshape(Hess, n, n);
    Hess = Hess + Hess.'; % ((M.*A)*resid) part is ready
    
    Hess = Hess + (grad_resid')*grad_resid; % Hess is ready
    
    % new iteration    
    du = Hess \ Grad;
%     du = pcg(Hess, Grad,tol, 100);
    du_norm = norm(du);
    u = u-du;
    
    fprintf('iter: %d, err_old: %3.3e, du: %3.3e\n', i, err_old, du_norm);
    if (err_old<tol)
        break;
    end;
end;


end