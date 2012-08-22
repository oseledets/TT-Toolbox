function [y]=fort_mvk4(A,x,eps,varargin)

d = x.d;
nswp = 10;
kickrank = 5;
rmax = 256;
y0 = [];
verb = 1;
for i=1:2:length(varargin)-1
    switch lower(varargin{i})
        case 'nswp'
            nswp=varargin{i+1};
        case 'rmax'
            rmax=varargin{i+1};
        case 'y0'
            y0=varargin{i+1};
        case 'kickrank'
            kickrank=varargin{i+1};
        case 'verb'
            verb=varargin{i+1};            

        otherwise
            error('Unrecognized option: %s\n',varargin{i});
    end
end

if (isempty(y0))
    y0 = tt_rand(A.n, d, 2);
end;

[ry,cry]=fort_mvk4_mex(int64(d), int64(A.n), int64(A.m), ...
    int64(x.r), int64(A.r), A.core, x.core, y0.core, int64(y0.r), eps, ...
    int64(rmax), int64(nswp), int64(kickrank), int64(verb));
% fprintf('fort_mvk4 done.\n');

y = y0;
ry = double(ry);
y.r = ry;
y.ps = cumsum([1; ry(1:d).*A.n.*ry(2:d+1)]);
y.core = cry;

end
