function [X,Y]=cg_trunc(W, F, r, tol, maxit, X0, Y0)

m1 = size(F,1);
m2 = size(F,2);
n1 = size(W{1},1);
n2 = size(W{2},1);
R = size(W{1},3);
if (nargin<5)||(isempty(maxit))
    maxit=100;
end;
gamma = 0e-10;

% RHS = W*F
core1 = reshape(permute(W{1}, [3 1 2]), R*n1, m1);
rhs = core1*F; % size R*n1, m2
rhs = reshape(rhs, R, n1, m2);
rhs = reshape(permute(rhs, [3 1 2]), m2*R, n1);
core2 = reshape(W{2}, n2, m2*R);
rhs = core2*rhs; % size n2, n1
rhs = rhs.';

if (nargin<7)||(isempty(X0))||(isempty(Y0))
    [X,S,Y]=svd(F, 'econ');
    X=X(:,1:r);
    Y=conj(Y(:,1:r))*S(1:r,1:r);
%     X = [randn(m1,1), randn(m1,r-1)];
%     Y = [randn(m2,1), randn(m2,r-1)];
else
    X=X0;
    Y=Y0; 
    [X,rv]=qr(X,0);
    Y = Y*(rv.');
end;

X_good=X;
Y_good=Y;

for it=1:maxit           
    
    if (it==1)
        % Residual
        resid = reshape(bfun2(W, X, Y), n1, n2) - rhs;
%         core1 = reshape(permute(W{1}, [3 1 2]), R*n1, m1);
%         resid = core1*X; % size R*n1, r
%         resid = reshape(resid, R, n1*r);
%         core2 = reshape(W{2}, n2*m2, R);
%         resid = core2*resid; % size n2*m2, n1*r
%         resid = reshape(resid, n2, m2, n1, r);
%         resid = reshape(permute(resid, [3 1 2 4]), n1*n2, m2*r);
%         resid = resid*reshape(Y, m2*r, 1); % size n1*n2,1
%         resid = reshape(resid, n1, n2) - rhs;
        
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
    
    Grad_full = [Xgrad; Ygrad];
    
    % direction update
    if (it>1)
        beta = (Grad_full')*(Grad_full) / ((Grad_prev')*Grad_prev);
        dir = Grad_full + beta*dir_prev;
    else
        dir = Grad_full;
    end;
    
    Grad_prev = Grad_full;
    dir_prev = dir;
    
    Xdir = reshape(dir(1:m1*r), m1, r);
    Ydir = reshape(dir(m1*r+1:(m1+m2)*r), m2, r);
    
    % Compute polynomial for alpha
    WdXdY = bfun2(W, Xdir, Ydir);
    WXdY = bfun2(W, X, Ydir);
    WdXY = bfun2(W, Xdir, Y);
    coeffs = zeros(1,4);
    coeffs(1) = dot(WdXdY, WdXdY)*2; % alpha^3
    coeffs(2) = dot(WdXdY, WXdY+WdXY)*3; % alpha^2
    coeffs(3) = dot(WdXdY, reshape(resid, n1*n2,1))*2 + dot(WXdY+WdXY,WXdY+WdXY); % alpha^1
    coeffs(4) = dot(WXdY+WdXY, reshape(resid, n1*n2,1)); % alpha^0
    % Frobenius matrix - for roots
    Frob_mtx = zeros(3,3);
    Frob_mtx(2,1)=1;
    Frob_mtx(3,2)=1;
    Frob_mtx(1,3)=-coeffs(4)/coeffs(1);
    Frob_mtx(2,3)=-coeffs(3)/coeffs(1);
    Frob_mtx(3,3)=-coeffs(2)/coeffs(1);
    alpha = eig(Frob_mtx);
    % Find actual minimum
%     errs = zeros(1,3);
%     for i=1:3
%         newX = X+alpha(i)*Xdir;
%         newY = Y+alpha(i)*Ydir;
%         resid = reshape(bfun2(W, newX, newY), n1, n2) - rhs;
%         errs(i) = norm(resid, 'fro')/norm(rhs, 'fro');
%     end;
% %     keyboard;
%     [~,i]=min(errs);
    i=1;
    alpha = real(alpha(i));
    
    X = X+alpha*Xdir;
    Y = Y+alpha*Ydir;
    
    resid = reshape(bfun2(W, X, Y), n1, n2) - rhs;
    err = norm(resid, 'fro')/norm(rhs, 'fro');
    conv = err_old/err;
    if (err<=err_min)
        X_good = X;
        Y_good = Y;
        err_min = err;
    end;
    err_old = err;    
    
    fprintf('iter %d: err = %3.3e, conv.fact = %3.5f\n', it, err, conv);
    
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