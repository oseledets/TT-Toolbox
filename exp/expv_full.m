function [u]=expv_full(A, u0, tol, varargin)
% function [u]=expv_full(A, u0, tol, [m, rescale, verb])

m = 20;
rescale = true;
verb = 0;
for i=1:2:length(varargin)-1
    switch lower(varargin{i})
        case 'm'
            m=varargin{i+1};
        case 'rescale'
            rescale=varargin{i+1};
        case 'verb'
            verb=varargin{i+1};
        otherwise
            error('expv_full: unknown option %s\n', varargin{i});
    end;
end;

atype = 0;
if (isa(A, 'function_handle'))
    atype = 1;
end;

n = size(u0,1);

if (rescale)
    v = randn(n, 1);
    v = v/norm(v);
    for i=1:5
        v = mvfun(A,atype,v);
        nrmA = norm(v);
        v = v/nrmA;
    end;
    nrmA = ceil(nrmA);
    if (verb>0)
        fprintf('exp-full: found |A|=%3.3e\n', nrmA);
    end;
else
    nrmA = 1;
end;

for t=1:nrmA    
    u = u0;
    nrm0 = norm(u0);
    du = u0;
    for i=2:m
        du = mvfun(A,atype,du);
        du = du/((i-1)*nrmA);
        u = u+du;
        if (verb>0)
            fprintf('exp-full: t=%d, i=%d, |du|=%3.3e\n', t, i, norm(du)/nrm0);
        end;
        if (norm(du)/nrm0<tol)
            break;
        end;
    end;
    u0 = u;
end;

end

function [y]=mvfun(A,atype,x)
if (atype==1)
    y = A(x);
else
    y = A*x;
end;
end


