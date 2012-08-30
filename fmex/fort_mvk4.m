function [y]=fort_mvk4(A,x,eps,varargin)

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

rx = x.r;
d = x.d;
tailranks = [false,false];
if (rx(1)~=1)
    e1 = tt_tensor;
    e1.d=1;
    e1.n = [rx(1)];
    e1.r = [1; rx(1)];
    e1.ps = [1; rx(1)^2+1];
    e1.core = reshape(eye(rx(1)), rx(1)^2, 1);
    x = kron(e1, x);
    
    A = kron(tt_matrix(eye(rx(1))), A);
    tailranks(1)=true;
end;
if (rx(d+1)~=1)
    e1 = tt_tensor;
    e1.d=1;
    e1.n = [rx(d+1)];
    e1.r = [rx(d+1); 1];
    e1.ps = [1; rx(d+1)^2+1];
    e1.core = reshape(eye(rx(d+1)), rx(d+1)^2, 1);
    x = kron(x, e1);
    
    A = kron(A, tt_matrix(eye(rx(d+1))));
    tailranks(2)=true;
end;

d = x.d;
n = A.n;

if (isempty(y0))
    y0 = tt_rand(n, d, 2);
end;

[ry,cry]=fort_mvk4_mex(int64(d), int64(n), int64(A.m), ...
    int64(x.r), int64(A.r), A.core, x.core, y0.core, int64(y0.r), eps, ...
    int64(rmax), int64(nswp), int64(kickrank), int64(verb));
% fprintf('fort_mvk4 done.\n');

ry = double(ry);
ps = cumsum([1; ry(1:d).*n.*ry(2:d+1)]);
if (tailranks(1))
    cr1 = reshape(cry(1:ps(2)-1), n(1), ry(2));
    cr2 = reshape(cry(ps(2):ps(3)-1), ry(2), n(2)*ry(3));
    cr2 = cr1*cr2;
    cry = [cr2(:); cry(ps(3):ps(d+1)-1)];
    ry = ry(2:d+1);
    ry(1) = n(1);
    n = n(2:d);
    d = d-1;
    ps = cumsum([1; ry(1:d).*n.*ry(2:d+1)]);
end;
if (tailranks(2))
    cr1 = reshape(cry(ps(d-1):ps(d)-1), ry(d-1)*n(d-1), ry(d));
    cr2 = reshape(cry(ps(d):ps(d+1)-1), ry(d), n(d));
    cr2 = cr1*cr2;
    cry = [cry(1:ps(d-1)-1); cr2(:)];
    d = d-1;
    ry = ry(1:d+1);
    ry(d+1) = n(d+1);
    n = n(1:d);
    ps = cumsum([1; ry(1:d).*n.*ry(2:d+1)]);
end;

y = y0;
y.r = ry;
y.core = cry;
y.ps = ps;
y.n = n;
y.d = d;

end
