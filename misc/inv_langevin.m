function [t]=inv_langevin(l, varargin)
% function [t]=inv_langevin(l, varargin)
% Computes the Inverse Langevin function in the nodes l via the Newton
% iterations.
% l should be in (-1,1), varargin are:
%   maxit - maximal number of Newton iterations, default 50
%   tol - stopping tolerance for ||L(t)-l||_2, default 1e-13
%   t0 - initial guess, default 0.5
%   verb - verbosity level:
%       0 - no messages (default)
%       1 - final information
%       2 - iteration information

tol = 1e-13;
maxit = 50;
t = 0.5*ones(size(l));
verb=0;
for i=1:2:length(varargin)-1
    switch lower(varargin{i})
        case 'maxit'
            maxit=varargin{i+1};
        case 'tol'
            tol=varargin{i+1};
        case 't0'
            t=varargin{i+1};            
        case 'verb'
            verb=varargin{i+1};                        
        otherwise
            error('Unrecognized option: %s\n',varargin{i});
    end
end


for i=1:maxit
    t = t-(coth(t)-1./t-l)./(-(sinh(t)).^(-2)+t.^(-2));
    t(l==0)=0;
    resid = coth(t)-1./t-l;
    resid(l==0)=0;
    normF = norm(resid, 'fro');
    if (verb>1)
        fprintf('iters: %d, |F|=%3.3e\n', i, normF);
    end;
    if (normF<tol)
        break;
    end;
end;

if (verb>0)
    fprintf('iters: %d, |F|=%3.3e\n', i, normF);
end;

end
