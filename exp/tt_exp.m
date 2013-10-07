function [y]=tt_exp(x, eps, varargin)
%Computation of the pointwise exponential in TT format
%   [Y]=TT_EXP(X,EPS,VARARGIN) This function computes pointwise exponential
%   using scaling and squaring method of the TT-vector X. EPS is the
%   accuracy parameter, varargins: 
%       N is the number of summand for the local Taylor
%           series (N=10 by default, usually enough). 
%       RMAX is the TT-rank bound.


% nrm = norm(x);
if (isa(x, 'tt_tensor'))
    nrm = tt_max_abs(x);
else
    nrm = tt_max_abs(qtttucker_to_linqtt(x, eps));
end;
n0 = floor(max(log2(nrm), 0))+1;
x = x./(2^n0);
if (isa(x, 'tt_tensor'))
    ons = tt_ones(x.n);
else
    ons = [];
    for i=1:(x.core.d)
        curons = tt_ones(x.tuck{i}.n);
        curons = qtt_tucker(curons, x.tuck{i}.d, eps);
        ons = kron(ons, curons);
    end;
end;
y = ons;

N = 10;
hdm = 'svd';
rmax = Inf;
epst = eps;
while (length(varargin)==1)
    varargin = varargin{1};
end;
for i=1:2:length(varargin)-1
    switch lower(varargin{i})
        case 'n'
            N=varargin{i+1};
        case 'rmax'
            rmax=varargin{i+1};
        case 'epst'
            epst=varargin{i+1};
        case 'hdm'
            hdm=varargin{i+1};
            
        otherwise
            error('Unrecognized option: %s\n',varargin{i});
    end
end

for k=(N-1):-1:1
% for k=1:N-1
    y=ons+(y.*x)/k;
    y=round(y,epst,rmax);
end

for k=1:n0
    if (strcmp(hdm, 'svd'))
        y=round(y.*y,eps*(0.5^(n0-k)),rmax);
        fprintf('squaring %d\n', k);
    else
        if (isa(x, 'tt_tensor'))
%             y = mvk3(diag(y), y, eps, 'nswp', 20, 'kickrank', 2);
            y = tt_mvk4(diag(y), y, eps, 'nswp', 20);
        else
%             y = mvrk(diag(y), y, eps, 'nswp', 20, 'kickrank', 2);
            y = mvrk2(diag(y), y, eps, 'nswp', 20);
        end;
    end;
end

end