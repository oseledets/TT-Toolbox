function [u,s,r]=svdgram(A,tol)
% Highly experimental acceleration of SVD/QR using Gram matrix.
% Use with caution for m>>n only!
%   function [u,s,r]=svdgram(A,[tol])
%   u is the left singular factor of A, 
%   s is the singular values (vector!), 
%   r has the meaning of diag(s)*v'.
%   if tol is given, performs the truncation with Fro-threshold.

R2 = A'*A;
[u,s,v]=svd(R2, 'econ');
u = A*v;
s = sum(u.^2, 1);
s = sqrt(s.');
if (nargin>1)&&(~isempty(tol))
    p = my_chop2(s, norm(s)*tol);
    u = u(:,1:p);
    s = s(1:p);
    v = v(:,1:p);
end;
u = u*spdiags(1./s, 0, numel(s), numel(s)); 
if (issparse(u))
    u = full(u);
end;
r = diag(s)*v';

% Run chol for reortogonalization.
% It will stop if the matrix will be singular.
% Fortunately, it means rank truncation with eps in our business.
if (s(1)/s(end)>1e7)
    p = 1;
    while (p>0)
        R2 = u'*u;
        [R,p] = chol(R2);
        if (p>0)
            u = u(:,1:p-1);
            s = s(1:p-1);
            r = r(1:p-1,:);
        end;
        iR = inv(R);
        u = u*iR;
        r = R*r;
    end;
end;

end