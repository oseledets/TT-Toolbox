function [y,z]=amen_mv(A, x, tol, varargin)
%Approximate the matrix-by-vector via the AMEn iteration
%   [y,z]=amen_mv(A, x, tol, varargin)
%   Attempts to approximate the y = A*x
%   with accuracy TOL using the AMEn+ALS iteration.
%   Matrix A has to be given in the TT-format, right-hand side x should be
%   given in the TT-format also. 
%
%   Options are provided in form
%   'PropertyName1',PropertyValue1,'PropertyName2',PropertyValue2 and so
%   on. The parameters are set to default (in brackets in the following)
%   The list of option names and default values are:
%       o y0 - initial approximation to Ax [rand rank-2]
%       o nswp - maximal number of sweeps [20]
%       o verb - verbosity level, 0-silent, 1-sweep info, 2-block info [1]
%       o kickrank - compression rank of the error, 
%         i.e. enrichment size [3]
%       o init_qr - perform QR of the input (save some time in ts, etc) [true]
%       o renorm - Orthog. and truncation methods: direct (svd,qr) or gram
%         (apply svd to the gram matrix, faster for m>>n) [direct]
%       o fkick - Perform solution enrichment during forward sweeps [false]
%         (rather questionable yet; false makes error higher, but "better
%         structured": it does not explode in e.g. subsequent matvecs)
%       o z0 - initial approximation to the error Ax-y [rand rank-kickrank]
%
%
%********
%   For description of adaptive ALS please see
%   Sergey V. Dolgov, Dmitry V. Savostyanov,
%   Alternating minimal energy methods for linear systems in higher dimensions. 
%   Part I: SPD systems, http://arxiv.org/abs/1301.6068,
%   Part II: Faster algorithm and application to nonsymmetric systems, http://arxiv.org/abs/1304.1222
%
%   Use {sergey.v.dolgov, dmitry.savostyanov}@gmail.com for feedback
%********
%
% TT-Toolbox 2.2, 2009-2012
%
%This is TT Toolbox, written by Ivan Oseledets et al.
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------


if (~isempty(varargin))
    v1 = varargin{1};
    if (isa(v1, 'cell'))
        varargin=v1;
    end;
end;

nswp = 20;
kickrank = 4;
kickrank2 = 0;
verb = 1;
y = [];
z = [];
init_qr = true;
renorm = 'direct';
% renorm = 'gram';
fkick = false;
% fkick = true;

for i=1:2:length(varargin)-1
    switch lower(varargin{i})
        case 'nswp'
            nswp=varargin{i+1};
        case 'y0'
            y=varargin{i+1};
        case 'z0'
            z=varargin{i+1};            
        case 'verb'
            verb=varargin{i+1};
        case 'kickrank'
            kickrank=varargin{i+1};
        case 'kickrank2'
            kickrank2=varargin{i+1};             
        case 'init_qr'
            init_qr = varargin{i+1};
        case 'renorm'
            renorm = varargin{i+1};
        case 'fkick'
            fkick = varargin{i+1};            
            
        otherwise
            warning('Unrecognized option: %s\n',varargin{i});
    end
end

if (isa(x, 'tt_tensor'))
    d = x.d;
    m = x.n;
    rx = x.r;
    x = core2cell(x);
    vectype = 1; % tt_tensor
else
    d = numel(x);
    m = zeros(d,1);
    rx = ones(d+1,1);
    for i=1:d
        m(i) = size(x{i}, 2);
        rx(i+1) = size(x{i}, 3);
    end;
    vectype = 0; % cell
end;

if (isa(A, 'tt_matrix'))
    n = A.n;
    ra = A.r;
    A = core2cell(A);
    % prepare A for fast ALS-mv
    for i=1:d
        A{i} = reshape(A{i}, ra(i)*n(i), m(i)*ra(i+1));
    end;
    atype = 1; % tt_matrix
% Alternative: A is a cell of cell: sparse canonical format
elseif (isa(A, 'cell'))
    n = zeros(d, 1);
    for i=1:d
        n(i) = size(A{i}{1}, 1);
    end;
    ra = numel(A{1});
    atype = 0; % cell
end;

if (isempty(y))
    y = tt_rand(n, d, 2);
    y = core2cell(y);
else
    if (isa(y, 'tt_tensor'))
        y = core2cell(y);
    end;
end;
ry = ones(d+1,1);
for i=1:d
    ry(i+1) = size(y{i}, 3);
end;

if (kickrank+kickrank2>0)
    if (isempty(z))
        z = tt_rand(n,d,kickrank+kickrank2);
        rz = z.r;
        z = core2cell(z);
    else
        if (isa(z, 'tt_tensor'))
            z = core2cell(z);
        end;
        rz = ones(d+1,1);
        for i=1:d
            rz(i+1) = size(z{i}, 3);
        end;
    end;
    
    phizax = cell(d+1,1); 
    if (atype==1)
        phizax{1}=1; phizax{d+1}=1;
    else
        phizax{1}=num2cell(ones(1,ra)); phizax{d+1}=num2cell(ones(1,ra));
    end;
    phizy = cell(d+1,1); phizy{1}=1; phizy{d+1}=1;
end;

phiyax = cell(d+1,1);
if (atype==1)
    phiyax{1}=1; phiyax{d+1}=1;
else
    phiyax{1}=num2cell(ones(1,ra)); phiyax{d+1}=num2cell(ones(1,ra));
end;

nrms = ones(d,1);

% Initial ort
for i=1:d-1
    if (init_qr)
        cr = reshape(y{i}, ry(i)*n(i), ry(i+1));
        if (strcmp(renorm, 'gram'))&&(ry(i)*n(i)>5*ry(i+1))
            [cr,s,R]=svdgram(cr);
        else
            [cr,R]=qr(cr, 0);
        end;
        nrmr = norm(R, 'fro');
        if (nrmr>0)
            R = R/nrmr;
        end;        
        cr2 = reshape(y{i+1}, ry(i+1), n(i+1)*ry(i+2));
        cr2 = R*cr2;
        ry(i+1) = size(cr, 2);
        y{i} = reshape(cr, ry(i), n(i), ry(i+1));
        y{i+1} = reshape(cr2, ry(i+1), n(i+1), ry(i+2));
    end;
    [phiyax{i+1},nrms(i)] = compute_next_Phi(phiyax{i}, y{i}, A{i}, x{i}, 'lr');
    
    if (kickrank+kickrank2>0)
        cr = reshape(z{i}, rz(i)*n(i), rz(i+1));
        if (strcmp(renorm, 'gram'))&&(rz(i)*n(i)>5*rz(i+1))
            [cr,s,R]=svdgram(cr);
        else
            [cr,R]=qr(cr, 0);
        end;
        nrmr = norm(R, 'fro');
        if (nrmr>0)
            R = R/nrmr;
        end;        
        cr2 = reshape(z{i+1}, rz(i+1), n(i+1)*rz(i+2));
        cr2 = R*cr2;
        rz(i+1) = size(cr, 2);
        z{i} = reshape(cr, rz(i), n(i), rz(i+1));
        z{i+1} = reshape(cr2, rz(i+1), n(i+1), rz(i+2));
        phizax{i+1} = compute_next_Phi(phizax{i}, z{i}, A{i}, x{i}, 'lr', nrms(i));
        phizy{i+1} = compute_next_Phi(phizy{i}, z{i}, [], y{i}, 'lr');
    end;
end;

i = d;
dir = -1;
swp = 1;
max_dx = 0;

while (swp<=nswp)
    % Project the MatVec generating vector
    crx = reshape(x{i}, rx(i)*m(i)*rx(i+1), 1);
    cry = bfun3(phiyax{i}, A{i}, phiyax{i+1}, crx);
    nrms(i) = norm(cry, 'fro');
    % The main goal is to keep y{i} of norm 1
    if (nrms(i)>0)
        cry = cry/nrms(i);
    else
        nrms(i)=1;
    end;    
    y{i} = reshape(y{i}, ry(i)*n(i)*ry(i+1), 1);
    dx = norm(cry-y{i});
    max_dx = max(max_dx, dx);
    
    % Truncation and enrichment
    if ((dir>0)&&(i<d))
        cry = reshape(cry, ry(i)*n(i), ry(i+1));
        if (strcmp(renorm, 'gram'))
            [u,s,v]=svdgram(cry, tol/sqrt(d));
            v = v.';
            r = size(u,2);
        else
            [u,s,v]=svd(cry, 'econ');
            s = diag(s);
            r = my_chop2(s, tol*norm(s)/sqrt(d));
            u = u(:,1:r);
            v = conj(v(:,1:r))*diag(s(1:r));
        end;        
        
        % Prepare enrichment, if needed
        if (kickrank+kickrank2>0)
            cry = u*v.';
            cry = reshape(cry, ry(i)*n(i), ry(i+1));
            % For updating z
            crz = bfun3(phizax{i}, A{i}, phizax{i+1}, crx);
            crz = reshape(crz, rz(i)*n(i), rz(i+1));            
            ys = cry*phizy{i+1};
            yz = reshape(ys, ry(i), n(i)*rz(i+1));
            yz = phizy{i}*yz;
            yz = reshape(yz, rz(i)*n(i), rz(i+1));
            crz = crz/nrms(i) - yz;            
            nrmz = norm(crz,'fro');
            if (kickrank2>0)
                [crz,~,~]=svd(crz, 'econ');
                crz = crz(:, 1:min(size(crz,2), kickrank));
                crz = [crz, randn(rz(i)*n(i), kickrank2)];
            end;
            % For adding into solution
            if (fkick)
                crs = bfun3(phiyax{i}, A{i}, phizax{i+1}, crx);
                crs = reshape(crs, ry(i)*n(i), rz(i+1));
                crs = crs/nrms(i) - ys;
                u = [u,crs];
                if (strcmp(renorm, 'gram'))&&(ry(i)*n(i)>5*(ry(i+1)+rz(i+1)))
                    [u,s,R]=svdgram(u);
                else
                    [u,R]=qr(u, 0);
                end;
                v = [v, zeros(ry(i+1), rz(i+1))];
                v = v*R.';
                r = size(u, 2);
            end;
        end;
        y{i} = reshape(u, ry(i), n(i), r);
        
        cr2 = reshape(y{i+1}, ry(i+1), n(i+1)*ry(i+2));
        v = reshape(v, ry(i+1), r);
        cr2 = v.'*cr2;
        y{i+1} = reshape(cr2, r, n(i+1), ry(i+2));
        
        ry(i+1) = r;
        
        [phiyax{i+1}, nrms(i)] = compute_next_Phi(phiyax{i}, y{i}, A{i}, x{i}, 'lr');
        
        if (kickrank+kickrank2>0)
            if (strcmp(renorm, 'gram'))&&(rz(i)*n(i)>5*rz(i+1))
                [crz,s,R]=svdgram(crz);
            else
                [crz,R]=qr(crz, 0);
            end;
            rz(i+1) = size(crz, 2);
            z{i} = reshape(crz, rz(i), n(i), rz(i+1));
            % z{i+1} will be recomputed from scratch in the next step
            
            phizax{i+1} = compute_next_Phi(phizax{i}, z{i}, A{i}, x{i}, 'lr', nrms(i));
            phizy{i+1} = compute_next_Phi(phizy{i}, z{i}, [], y{i}, 'lr');
        end;
    elseif ((dir<0)&&(i>1))
        cry = reshape(cry, ry(i), n(i)*ry(i+1));
        if (strcmp(renorm, 'gram'))
            [v,s,u]=svdgram(cry.', tol/sqrt(d));
            u = u.';
            r = size(v,2);
        else
            [u,s,v]=svd(cry, 'econ');
            s = diag(s);
            r = my_chop2(s, tol*norm(s)/sqrt(d));
            v = conj(v(:,1:r));
            u = u(:,1:r)*diag(s(1:r));
        end;
        
        % Prepare enrichment, if needed
        if (kickrank+kickrank2>0)
            cry = u*v.';
            cry = reshape(cry, ry(i), n(i)*ry(i+1));
            % For updating z
            crz = bfun3(phizax{i}, A{i}, phizax{i+1}, crx);
            crz = reshape(crz, rz(i), n(i)*rz(i+1));            
            ys = phizy{i}*cry;
            yz = reshape(ys, rz(i)*n(i), ry(i+1));
            yz = yz*phizy{i+1};
            yz = reshape(yz, rz(i), n(i)*rz(i+1));
            crz = crz/nrms(i) - yz;
            nrmz = norm(crz,'fro');
            if (kickrank2>0)
                [~,~,crz]=svd(crz, 'econ');
                crz = crz(:, 1:min(size(crz,2), kickrank))';
                crz = [crz; randn(kickrank2, n(i)*rz(i+1))];
            end;            
            % For adding into solution
            crs = bfun3(phizax{i}, A{i}, phiyax{i+1}, crx);
            crs = reshape(crs, rz(i), n(i)*ry(i+1));
            crs = crs/nrms(i) - ys;
            v = [v,crs.'];
            if (strcmp(renorm, 'gram'))&&(n(i)*ry(i+1)>5*(ry(i)+rz(i)))
                [v,s,R]=svdgram(v);
            else
                [v,R]=qr(v, 0);
            end;
            u = [u, zeros(ry(i), rz(i))];
            u = u*R.';
            r = size(v, 2);
        end;
        cr2 = reshape(y{i-1}, ry(i-1)*n(i-1), ry(i));
        cr2 = cr2*u;
        y{i-1} = reshape(cr2, ry(i-1), n(i-1), r);
        y{i} = reshape(v.', r, n(i), ry(i+1));
        
        ry(i) = r;
        
        [phiyax{i},nrms(i)] = compute_next_Phi(phiyax{i+1}, y{i}, A{i}, x{i}, 'rl');
        
        if (kickrank+kickrank2>0)
            if (strcmp(renorm, 'gram'))&&(n(i)*rz(i+1)>5*rz(i))
                [crz,s,R]=svdgram(crz.');
            else
                [crz,R]=qr(crz.', 0);
            end;
            rz(i) = size(crz, 2);
            z{i} = reshape(crz.', rz(i), n(i), rz(i+1));
            % don't update z{i-1}, it will be recomputed from scratch
            
            phizax{i} = compute_next_Phi(phizax{i+1}, z{i}, A{i}, x{i}, 'rl', nrms(i));
            phizy{i} = compute_next_Phi(phizy{i+1}, z{i}, [], y{i}, 'rl');
        end;
    end;
    
    if (verb>1)
        fprintf('amen-mv: swp=[%d,%d], dx=%3.3e, r=%d, |y|=%3.3e, |z|=%3.3e\n', swp, i, dx, r, norm(cry,'fro'), nrmz);
    end;
    
    % Stopping or reversing
    if ((dir>0)&&(i==d))||((dir<0)&&(i==1))
        if (verb>0)
            fprintf('amen-mv: swp=%d{%d}, max_dx=%3.3e, max_r=%d\n', swp, (1-dir)/2, max_dx, max(ry));
        end;
        if ((max_dx<tol)||(swp==nswp))&&(dir>0)
            break;
        else
            % We are at the terminal block
            y{i} = reshape(cry, ry(i), n(i), ry(i+1));
            if (dir>0); swp = swp+1; end;
        end;
        max_dx = 0;
        dir = -dir;
    else
        i = i+dir;
    end;
end;
% if (dir>0)
    y{d} = reshape(cry, ry(d), n(d), ry(d+1));
% else
%     y{1} = reshape(cry, ry(1), n(1), ry(2));
% end;

% Distribute norms equally...
nrms = exp(sum(log(nrms))/d);
% ... and plug them into y
for i=1:d
    y{i} = y{i}*nrms;
end;


if (vectype==1)
    y = cell2core(tt_tensor, y);
    z = cell2core(tt_tensor, z);
end;
end


% new
function [Phi,nrm] = compute_next_Phi(Phi_prev, x, A, y, direction,extnrm)
% Performs the recurrent Phi (or Psi) matrix computation
% Phi = Phi_prev * (x'Ay).
% If direction is 'lr', computes Psi
% if direction is 'rl', computes Phi
% A can be empty, then only x'y is computed.

% Phi1: rx1, ry1, ra1, or {rx1, ry1}_ra, or rx1, ry1
% Phi2: ry2, ra2, rx2, or {ry2, rx2}_ra, or ry2, rx2

if (nargin<6)
    extnrm = [];
end;


rx1 = size(x,1); n = size(x,2); rx2 = size(x,3);
ry1 = size(y,1); m = size(y,2); ry2 = size(y,3);
if (~isempty(A))
    if (isa(A, 'cell'))
        % A is a canonical block
        ra = numel(A); 
    else
        % Just full format
        ra1 = size(A,1)/n; ra2 = size(A,2)/m;
    end;
else
    ra1 = 1; ra2 = 1;
end;

if (isa(Phi_prev, 'cell'))
    Phi = cell(ra, 1);
    if (nargout>1)
        nrm = 0;
    end;
    if (strcmp(direction, 'lr'))
        %lr: Phi1
        x = reshape(x, rx1, n*rx2);
        y = reshape(y, ry1*m, ry2);
        for i=1:ra
            Phi{i} = x'*Phi_prev{i};
            Phi{i} = reshape(Phi{i}, n, rx2*ry1);            
            Phi{i} = Phi{i}.';            
            Phi{i} = Phi{i}*A{i};
            Phi{i} = reshape(Phi{i}, rx2, ry1*m);
            Phi{i} = Phi{i}*y;
            if (nargout>1)
                nrm = max(nrm, norm(Phi{i}, 'fro'));
            end;            
        end;
    else
        %rl: Phi2
        y = reshape(y, ry1, m*ry2);
        x = reshape(x, rx1*n, rx2);
        for i=1:ra
            Phi{i} = Phi_prev{i}*x';
            Phi{i} = reshape(Phi{i}, ry2*rx1, n);
            Phi{i} = Phi{i}*A{i};
            Phi{i} = Phi{i}.';
            Phi{i} = reshape(Phi{i}, m*ry2, rx1);
            Phi{i} = y*Phi{i};
            if (nargout>1)
                nrm = max(nrm, norm(Phi{i}, 'fro'));
            end;
        end;
    end;
    if (nargout>1)
        % Extract the scale to prevent overload
        if (nrm>0)
            for i=1:ra
                Phi{i} = Phi{i}/nrm;
            end;
        else
            nrm=1;
        end;
    elseif (~isempty(extnrm))
        % Override the normalization
        for i=1:ra
            Phi{i} = Phi{i}/extnrm;
        end;
    end;
else
    if (strcmp(direction, 'lr'))
        %lr: Phi1
        x = reshape(x, rx1, n*rx2);
        Phi = reshape(Phi_prev, rx1, ry1*ra1);
        Phi = x'*Phi;
        if (~isempty(A))
            Phi = reshape(Phi, n*rx2*ry1, ra1);
            Phi = Phi.';
            Phi = reshape(Phi, ra1*n, rx2*ry1);
            Phi = A.'*Phi;
            Phi = reshape(Phi, m, ra2*rx2*ry1);
        else
            Phi = reshape(Phi, n, rx2*ry1);
        end;
        Phi = Phi.';
        Phi = reshape(Phi, ra2*rx2, ry1*m);        
        y = reshape(y, ry1*m, ry2);
        Phi = Phi*y;
        if (~isempty(A))
            Phi = reshape(Phi, ra2, rx2*ry2);
            Phi = Phi.';
        end;
        Phi = reshape(Phi, rx2, ry2, ra2);
    else
        %rl: Phi2
        y = reshape(y, ry1*m, ry2);
        Phi = reshape(Phi_prev, ry2, ra2*rx2);
        Phi = y*Phi;
        if (~isempty(A))
            Phi = reshape(Phi, ry1, m*ra2*rx2);
            Phi = Phi.';
            Phi = reshape(Phi, m*ra2, rx2*ry1);
            Phi = A*Phi;
            Phi = reshape(Phi, ra1*n*rx2, ry1);
            Phi = Phi.';
        end;        
        Phi = reshape(Phi, ry1*ra1, n*rx2);
        x = reshape(x, rx1, n*rx2);
        Phi = Phi*x';
        if (~isempty(A))
            Phi = reshape(Phi, ry1, ra1, rx1);
        else
            Phi = reshape(Phi, ry1, rx1);
        end;
    end;
    
    if (nargout>1)
        % Extract the scale to prevent overload
        nrm = norm(Phi(:), 'fro');
        if (nrm>0)
            Phi = Phi/nrm;
        else
            nrm=1;
        end;
    elseif (~isempty(extnrm))
        % Override the normalization by the external one
        Phi = Phi/extnrm;
    end;
end;

end


% new
function [y]=bfun3(Phi1, A, Phi2, x)
b = size(x,2);

if (isa(A, 'cell'))
    ra = numel(A);
    ry1 = size(Phi1{1},1);
    rx1 = size(Phi1{1},2);
    ry2 = size(Phi2{1},2);
    rx2 = size(Phi2{1},1);    
    n = size(A{1},1);
    m = size(A{1},2);
    
    y = zeros(ry1*n*ry2, b);
    for i=1:ra
        cy = reshape(x.', b*rx1*m, rx2);
        cy = cy*Phi2{i};
        cy = reshape(cy, b*rx1, m*ry2);
        cy = cy.';
        cy = reshape(cy, m, ry2*b*rx1);
        cy = A{i}*cy;
        cy = reshape(cy, n*ry2*b, rx1);
        cy = cy*Phi1{i}.';
        cy = cy.';
        cy = reshape(cy, ry1*n*ry2, b);
        y = y+cy;
    end;
else
    % Phi1: ry1, rx1, ra1
    ry1 = size(Phi1,1);
    rx1 = size(Phi1,2);
    ra1 = size(Phi1,3);
    % Phi2: rx2, ra2, ry2
    ry2 = size(Phi2,3);
    rx2 = size(Phi2,1);
    ra2 = size(Phi2,2);
    n = size(A,1)/ra1;
    m = size(A,2)/ra2;
    
    y = reshape(x.', b*rx1*m, rx2);
    Phi2 = reshape(Phi2, rx2, ra2*ry2);
    y = y*Phi2;
    y = reshape(y, b*rx1, m*ra2*ry2);
    y = y.';
    y = reshape(y, m*ra2, ry2*b*rx1);
    y = A*y;
    y = reshape(y, ra1*n*ry2*b, rx1);
    y = y.';
    y = reshape(y, rx1*ra1, n*ry2*b);
    Phi1 = reshape(Phi1, ry1, rx1*ra1);
    y = Phi1*y;
    y = reshape(y, ry1*n*ry2, b);
end;
end
