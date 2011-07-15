function [u]=cg_solve(B, M, A, f, tol, u0)

maxit=400;
n = size(B,1);
if (nargin<7)||(isempty(u0))
    u0 = ones(n,1);
end;

u = u0;
resid = B*u + (M*u).*(A*u) - f;
err_old = norm(resid);
grad_old = zeros(size(f));
dir = zeros(size(f));
for i=1:maxit
%     Bu = B*u;
    Mu = M*u;
    Au = A*u;
    
    diag_Au = spdiags(Au, 0, n, n);
    diag_Mu = spdiags(Mu, 0, n, n);
    
    grad_resid = B + diag_Mu*A + diag_Au*M; % [i,\hat{i}] - matrix
    Grad = (grad_resid')*resid;
    
    beta = dot(grad_old, grad_old);
    if (beta~=0)
        beta = dot(Grad,Grad)/beta;
    end;
    
    dir = beta*dir + Grad;
    mix1 = B*dir + (M*u).*(A*dir)+(A*u).*(M*dir);
    mix2 = (M*dir).*(A*dir);
    c0 = dot(resid,resid);
    c1 = 2*dot(mix1,resid);
    c2 = 2*dot(mix2,resid)+dot(mix1,mix1);
    c3 = 2*dot(mix1,mix2);
    c4 = dot(mix2,mix2);
    
    % Frobenius matrix - for roots
    Frob_mtx = zeros(4,4);
    Frob_mtx(2,1)=1;
    Frob_mtx(3,2)=1;
    Frob_mtx(4,3)=1;
    Frob_mtx(1,4)=-c0/c4;
    Frob_mtx(2,4)=-c1/c4;
    Frob_mtx(3,4)=-c2/c4;
    Frob_mtx(4,4)=-c3/c4;
    alpha = eig(Frob_mtx);    
    alpha = real(alpha);
    [vv,j]=min(c4*alpha.^4 + c3*alpha.^3 + c2*alpha.^2 + c1*alpha + c0);
    alpha = alpha(j);
%     % A very tiny 1D newton for alpha
%     alpha = -1e-3;
%     for j=1:20
%         agrad = 4*c4*alpha^3 + 3*c3*alpha^2 + 2*c2*alpha + c1;
%         ahess = 12*c4*alpha^2 + 6*c3*alpha + 2*c2;
%         alpha = alpha - agrad/ahess;
%         polyerr=c4*alpha^4 + c3*alpha^3 + c2*alpha^2 + c1*alpha + c0;
%         if (polyerr<1e-13)
%             break;
%         end;
%     end;
    
    u = u + alpha*dir;
    resid = B*u + (M*u).*(A*u) - f;
    err = norm(resid)/norm(f);
    conv = err_old/err;
    
%     fprintf('---cg solve: iter [%d,%d], err %3.3e, conv: %3.6f, poly_err: %3.3e\n',iout,i,err,conv,polyerr);
    grad_old = Grad;
    if (err<tol)
        break;
    end;
    if (conv-1<1e-7)
        break;
    end;
end;

fprintf('---cg solve: iter [%d], err %3.3e\n',i,err);

end