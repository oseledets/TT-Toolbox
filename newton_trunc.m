function [X,Y]=newton_trunc(W, F, r, tol, X0, Y0)

m1 = size(F,1);
m2 = size(F,2);
n1 = size(W{1},1);
n2 = size(W{2},1);
R = size(W{1},3);
maxit=10;
gamma = 0e-10;

% RHS = W*F
core1 = reshape(permute(W{1}, [3 1 2]), R*n1, m1);
rhs = core1*F; % size R*n1, m2
rhs = reshape(rhs, R, n1, m2);
rhs = reshape(permute(rhs, [3 1 2]), m2*R, n1);
core2 = reshape(W{2}, n2, m2*R);
rhs = core2*rhs; % size n2, n1
rhs = rhs.';

if (nargin<6)||(isempty(X0))||(isempty(Y0))
    [X,S,Y]=svd(F, 'econ');
    X=X(:,1:r);
    Y=conj(Y(:,1:r))*S(1:r,1:r);
%     X = 1e-2*[randn(m1,1), randn(m1,r-1)];
%     Y = [randn(m2,1), randn(m2,r-1)];
else
    X=X0;
    Y=Y0;
    [X,rv]=qr(X,0);
    Y = Y*(rv.');    
end;


for it=1:maxit           
    
    if (it==1)
        % Residual
        resid = bfun2(W, X, Y);
%         core1 = reshape(permute(W{1}, [3 1 2]), R*n1, m1);
%         resid = core1*X; % size R*n1, r
%         resid = reshape(resid, R, n1*r);
%         core2 = reshape(W{2}, n2*m2, R);
%         resid = core2*resid; % size n2*m2, n1*r
%         resid = reshape(resid, n2, m2, n1, r);
%         resid = reshape(permute(resid, [3 1 2 4]), n1*n2, m2*r);
%         resid = resid*reshape(Y, m2*r, 1); % size n1*n2,1
        resid = reshape(resid, n1, n2) - rhs;
        
        err_old = norm(resid, 'fro')/norm(rhs, 'fro');
        err_min = err_old;
    end;
    
    % X-gradient
    core1 = conj(reshape(permute(W{1}, [3 2 1]), R*m1, n1));
    Xgrad = core1*resid; % size R*m1, n2
    Xgrad = reshape(Xgrad, R, m1, n2);
    Xgrad = reshape(permute(Xgrad, [3 1 2]), n2*R, m1);
    core2 = conj(reshape(permute(W{2}, [2 1 3]), m2, n2*R));
    Xgrad = core2*Xgrad; % size m2, m1
    grad_save = Xgrad; % save for Ygrad and Hess12 calculation
    Xgrad = (Y')*Xgrad; % size r, m1
    Xgrad = reshape(Xgrad.', m1*r, 1) + gamma*reshape(X, m1*r, 1);
    
    % Y-gradient
    Ygrad = (X')*(grad_save.'); % size r, m2
    Ygrad = reshape(Ygrad.', m2*r, 1) + gamma*reshape(Y, m2*r, 1);
    
    % Hessian block 1,1 (X,X)
    core2 = reshape(permute(W{2}, [3 1 2]), R*n2, m2);
    Hess11 = core2*Y; % size R*n2, r
    core1 = reshape(W{1}, n1*m1, R);
    Hess11 = reshape(Hess11, R, n2*r);
    Hess11 = core1*Hess11; % size n1*m1, n2*r
    Hess11 = reshape(Hess11, n1, m1, n2, r);
    Hess11 = reshape(permute(Hess11, [1 3 2 4]), n1*n2, m1*r);
    Hess11 = (Hess11')*Hess11; % size m1*r, m1*r
    Hess11 = Hess11 + gamma*eye(m1*r);
    
    % Hessian block 1-2 (X,Y)
    core1 = reshape(permute(W{1}, [3 1 2]), R*n1, m1);
    Hess12 = core1*X; % size R*n1, r
    Hess12 = reshape(Hess12, R, n1, r);
    Hess12 = reshape(permute(Hess12, [2 1 3]), n1, R*r);
    core1 = conj(reshape(permute(W{1}, [3 2 1]), R*m1, n1));
    Hess12 = core1*Hess12; % size R*m1, R*r
    Hess12 = reshape(Hess12, R, m1, R, r);
    Hess12 = reshape(permute(Hess12, [3 4 1 2]), R, r*R*m1);
    core2 = reshape(W{2}, n2*m2, R);
    Hess12 = core2*Hess12; % size n2*m2, r*R*m1
    Hess12 = reshape(Hess12, n2, m2, r, R, m1);
    Hess12 = reshape(permute(Hess12, [1 4 2 3 5]), n2*R, m2*r*m1);
    core2 = conj(reshape(permute(W{2}, [2 1 3]), m2, n2*R));
    Hess12 = core2*Hess12; % size m2, m2*r*m1
    Hess12 = (Y')*Hess12; % size r, m2*r*m1
    Hess12 = reshape(Hess12, r, m2, r, m1);
    Hess12 = reshape(permute(Hess12, [4 1 2 3]), m1*r, m2*r);
    Ir = eye(r);
    Hess12 = Hess12 + kron(Ir, grad_save.');
    
    % Hessian block 2,2 (Y,Y)
    core1 = reshape(permute(W{1}, [3 1 2]), R*n1, m1);
    Hess22 = core1*X; % size R*n1, r
    core2 = reshape(W{2}, n2*m2, R);
    Hess22 = reshape(Hess22, R, n1*r);
    Hess22 = core2*Hess22; % size n2*m2, n1*r
    Hess22 = reshape(Hess22, n2, m2, n1, r);
    Hess22 = reshape(permute(Hess22, [3 1 2 4]), n1*n2, m2*r);
    Hess22 = (Hess22')*Hess22; % size m2*r, m2*r    
    Hess22 = Hess22 + gamma*eye(m2*r);
    
    % Full objects
    Hess = [Hess11, Hess12; Hess12', Hess22];
    Grad = [Xgrad; Ygrad];
    
    % Newton iteration
    sol = [reshape(X, m1*r, 1); reshape(Y, m2*r, 1)];

%     iHess = fgmres_selfprec(Hess, 1e-2, 1, 10);
%     dsol = Hess \ Grad;
%     dsol = gmres(@(v)(Hess*iHess*v), Grad, 50, 1e-10, 20);
%     dsol = iHess*dsol;
    dsol = gmres(@(v)(Hess*v), Grad, 50, 1e-10, 100);
%     dsol = pcg(Hess, Grad, 1e-10, 10000);
    sol = sol - dsol;
    
    X = reshape(sol(1:m1*r), m1, r);
    Y = reshape(sol((m1*r+1):((m1+m2)*r)), m2, r);
    
    % Residual
    resid = bfun2(W,X,Y);
%     core1 = reshape(permute(W{1}, [3 1 2]), R*n1, m1);
%     resid = core1*X; % size R*n1, r
%     resid = reshape(resid, R, n1*r);
%     core2 = reshape(W{2}, n2*m2, R);
%     resid = core2*resid; % size n2*m2, n1*r
%     resid = reshape(resid, n2, m2, n1, r);
%     resid = reshape(permute(resid, [3 1 2 4]), n1*n2, m2*r);
%     resid = resid*reshape(Y, m2*r, 1); % size n1*n2,1
    resid = reshape(resid, n1, n2) - rhs;
    
    err = norm(resid, 'fro')/norm(rhs, 'fro');
    conv = err_old/err;
    if (err<=err_min)
        X_good = X;
        Y_good = Y;
        err_min = err;
    end;
    err_old = err;    
    
    fprintf('iter %d: err = %3.3e, conv.fact = %3.5f, cond(Hess)=%3.3e\n', it, err, conv, cond(Hess));
    
    if (err<tol)
        break;
    end;
end;

X = X_good;
Y = Y_good;

end


function [Z]=bfun2(W, X, Y)
    m1 = size(X,1);
    m2 = size(Y,1);
    r = size(X,2);
    n1 = size(W{1},1);
    n2 = size(W{2},1);
    R = size(W{1},3);
    
    core1 = reshape(permute(W{1}, [3 1 2]), R*n1, m1);
    Z = core1*X; % size R*n1, r
    Z = reshape(Z, R*n1,r);
    Z = Z*(Y.'); % size R*n1, m2
    Z = reshape(Z, R, n1, m2);
    Z = reshape(permute(Z, [2 3 1]), n1, m2*R);
    core2 = reshape(W{2}, n2, m2*R);
    Z = Z*(core2.'); % size n1,n2
    Z = reshape(Z, n1*n2, 1);
end