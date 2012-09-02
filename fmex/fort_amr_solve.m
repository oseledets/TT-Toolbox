function [x]=fort_amr_solve(A,y,eps,varargin)

d = y.d;
nswp = 10;
kickrank = 4;
local_restart=40;
local_iters=2;
verb=1;
local_prec = 'n';
x0 = [];
for i=1:2:length(varargin)-1
    switch lower(varargin{i})
        case 'nswp'
            nswp=varargin{i+1};
        case 'x0'
            x0=varargin{i+1};
        case 'kickrank'
            kickrank=varargin{i+1};
        case 'local_restart'
            local_restart=varargin{i+1};
        case 'local_iters'
            local_iters=varargin{i+1};
        case 'local_prec'
            local_prec=lower(varargin{i+1});
        case 'verb'
            verb=varargin{i+1};

        otherwise
            error('Unrecognized option: %s\n',varargin{i});
    end
end

if (isempty(x0))
    x0 = tt_rand(A.n, d, 2);
end;

params = [nswp; kickrank; local_restart; local_iters; verb];
params = int64(params);

%  keyboard;

[rx,crx]=fort_amr_solve_mex(int64(d), int64(A.n), int64(A.m), ...
    int64(y.r), int64(A.r), A.core, y.core, x0.core, int64(x0.r), eps, params, local_prec);

x = x0;
rx = double(rx);
x.r = rx;
x.ps = cumsum([1; rx(1:d).*A.n.*rx(2:d+1)]);
x.core = crx;

end