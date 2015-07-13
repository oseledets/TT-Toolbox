function x_ = permute(x_, order, eps)
%Permute dimensions.
%   Y = permute(X, ORDER, EPS) permutes the dimensions of the TT-tensor X
%   according to ORDER, delivering a result at relative accuracy EPS. This
%   function is equivalent to 
%      Y = tt_tensor(permute(reshape(full(X), X.n'),ORDER), EPS) 
%   but avoids the conversion to the full format.
%
%
% Simon Etter, Summer 2015
% Seminar of Applied Mathematics, ETH Zurich
%
%
% TT-Toolbox 2.2, 2009-2012
%
%This is TT Toolbox, written by Ivan Oseledets et al. Institute of
%Numerical Mathematics, Moscow, Russia webpage:
%http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com 
%---------------------------

% This code basically performs insertion sort on the TT dimensions: 
%  for k = 2:d
%     Bubble the kth dimension to the right (according to ORDER) position in the first 1:k dimensions. 
%
% The current code could be optimised at the following places:
%  - Instead of initially orthogonalising with respect to the first two vertices, 
%    orthogonalise directly with respect to the first inversion. 
%  - When performing the SVD, check on which side of the current position the 
%    next swap will occur and directly establish the appropriate orthogonality
%    (current implementation always assumes we move left). 
% Both changes reduce the number of QRs by at most O(d) and are therefore likely
% to produce negligible speedup while rendering the code more complicated. 

% Parse input
x = core2cell(x_);
d = x_.d;
r = x_.r;
n = x_.n;
idx(order) = 1:length(x);
eps = eps/d^1.5; 
% ^Numerical evidence suggests that eps = eps/d may be sufficient, however I can only prove correctness
%  for this choice of global-to-local conversion factor.
if (length(order) < d) error('ORDER must have at least D elements for a D-dimensional tensor'); end;

% RL-orthogonalise x
for kk = d:-1:3
    [Q,R] = qr(reshape(x{kk}, [r(kk), n(kk)*r(kk+1)])',0);
    tr = min([r(kk), n(kk)*r(kk+1)]);
    x{kk} = reshape(Q', [tr, n(kk), r(kk+1)]);
    x{kk-1} = reshape(reshape(x{kk-1}, [r(kk-1)*n(kk-1),r(kk)])*R', [r(kk-1), n(kk-1), tr]);
    r(kk) = tr;
end
k = 1;
while (true)
    % Find next inversion
    nk = k;
    while (nk < d && idx(nk) < idx(nk+1)) nk = nk + 1; end;
    if (nk == d) break; end;

    % Move orthogonal centre there
    for kk = k:nk-1
        [Q,R] = qr(reshape(x{kk}, [r(kk)*n(kk), r(kk+1)]),0);
        tr = min([r(kk)*n(kk), r(kk+1)]);
        x{kk} = reshape(Q, [r(kk), n(kk), tr]);
        x{kk+1} = reshape(R*reshape(x{kk+1}, [r(kk+1), n(kk+1)*r(kk+2)]), [tr, n(kk+1), r(kk+2)]);
        r(kk+1) = tr;
    end
    k = nk;

    % Swap dimensions
    c = reshape(reshape(x{k}, [r(k)*n(k),r(k+1)])*reshape(x{k+1}, [r(k+1),n(k+1)*r(k+2)]), [r(k), n(k),n(k+1), r(k+2)]);
    c = permute(c, [1,3,2,4]);
    [U,S,V] = svd(reshape(c, [r(k)*n(k+1), n(k)*r(k+2)]), 'econ'); s = diag(S);
    r(k+1) = max(my_chop2(s, norm(s)*eps), 1);
    tmp = U*S; x{k} = reshape(tmp(:,1:r(k+1)), [r(k), n(k+1), r(k+1)]);
    x{k+1} = reshape(V(:,1:r(k+1))', [r(k+1), n(k), r(k+2)]);
    idx([k,k+1]) = idx([k+1,k]);
    n  ([k,k+1]) = n  ([k+1,k]);
    k = max(k-1, 1);
end

% Parse output
x_ = cell2core(x_,x);
end
