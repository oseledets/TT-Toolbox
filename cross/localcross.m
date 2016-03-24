function [u,v,I]=localcross(Y,tol)
% Full-pivoted cross for truncating one ket TT block instead of SVD

if (exist('localcross_mex','file')>0)&&(isreal(Y))
    % Check for a faster mex
    if (nargout>2)
        [u,v,I]=localcross_mex(Y,tol);
    else
        [u,v]=localcross_mex(Y,tol);
    end;
    return;
end;

[n,m,b]=size(Y);

minsz = min(n, m*b);
u = zeros(n, minsz);
v = zeros(minsz, m*b);
res = reshape(Y, n, m*b);
% Return also the indices
I = zeros(1, minsz);
val_max = max(abs(Y(:)));
for r=1:minsz
    res = reshape(res, [], 1);
    [val,piv]=max(abs(res));
    piv = tt_ind2sub([n, m*b], piv);
    if (val<=tol*val_max)
        break;
    end;
    res = reshape(res, n, m*b);
    u(:,r) = res(:, piv(2));
    v(r,:) = res(piv(1), :)/res(piv(1),piv(2));
    res = res - u(:,r)*v(r,:);
    I(r)=piv(1);
end;
% Return indices
r = find(I==0,1);
if (isempty(r))
    r=minsz;
else
    r=r-1;
end;
I = I(1:r);
u = u(:, 1:r);
v = v(1:r, :);
if (r==0)
    u = zeros(n,1);
    v = zeros(1,m);
    I = 1;
end;
% qr u, in case we don't have enrichment
[u,rv]=qr(u,0);
v = rv*v;
end

