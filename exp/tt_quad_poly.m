function [tt]=tt_quad_poly(B, f, tol)
% Generates the quadratic form (Lyapunov function).
% 
% function [TT]=tt_quad_poly(B, F, TOL)
% Generates the Tucker-core of V=x'Bx, and, optionally, the TT representation.
% B is the system matrix, optional parameters:
%   F - cell array of Tucker factors. Each of them should have the items
%   [1, x, x^2]. If F is omitted, returns only the Tucker core.
%   TOL is the SVD tolerance for the off-diagonal block decomposition
%
%   !!BYDLOCODE HERE!!: don't set large TOL!! TT-blocks are NOT orthogonal!!
%   Currently 1e-14 is set by default

if (nargin<3)||(isempty(tol))
    tol=1e-14;
end;

B = (B+B')*0.5; % The quadratic form uses only the symm part anyway

d = size(B,1);
tt = cell(d,1);

% First block
tt{1} = zeros(1, 3, 3);
tt{1}(1,3,1)=B(1,1);
tt{1}(1,2,2)=2;
tt{1}(1,1,3)=1;
% Previous U is 1, W is the first column of B
U_prev = 1;
W_prev = B(1,2:d);
r_prev = 1;

for i=2:d-1
    % Off-diagonal block decomp
    [U,si,wi]=svd(B(1:i, i+1:d), 'econ');
    if (norm(diag(si))==0)
        r=0;
    else
        r = my_chop2(diag(si), tol*norm(diag(si)));
    end;
    U = U(:, 1:r);
    W=si(1:r, 1:r)*wi(:,1:r)';
    
    Q = U_prev'*U(1:i-1,:);
    
    % Now, construct the k-th block
    tt{i} = zeros(r_prev+2, 3, r+2);
    tt{i}(1,1,1)=1;
    tt{i}(2:r_prev+1, 2, 1) = W_prev(:,1);
    tt{i}(r_prev+2, 3, 1) = B(i,i);
    tt{i}(2:r_prev+1, 1, 2:r+1) = reshape(Q, r_prev, 1, r);
    tt{i}(r_prev+2, 2, 2:r+1) = reshape(2*U(i, :), 1, 1, r);
    tt{i}(r_prev+2, 1, r+2) = 1;
    
    r_prev = r;
    U_prev = U;
    W_prev = W;
end;

% Now, the last block
tt{d} = zeros(r_prev+2, 3, 1);
tt{d}(1,1,1)=1;
tt{d}(2:r_prev+1,2,1) = W_prev;
tt{d}(r_prev+2,3,1) = B(d,d);

% OOP stuff
tt = cell2core(tt_tensor, tt);

if (nargin>=2)&&(~isempty(f))
    qt = qtt_tucker;
    qt.tuck = f;
    qt.core = tt;
    tt = qtttucker_to_tt(qt);
end;

end