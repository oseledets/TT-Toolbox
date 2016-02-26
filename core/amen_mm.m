function [X]=amen_mm(A, Y, tol, varargin)
%Approximate the matrix-by-matrix via the AMEn iteration
%   [X]=amen_mv(A, Y, tol, varargin)
%   Attempts to approximate the X = A*Y with accuracy TOL using the 
%   AMEn+ALS iteration. A is a n x m matrix, Y is a m x k matrix.
%   Matrices A,Y can be given in tt_matrix format, the output is tt_matrix.
%   Y can be tt_tensor, it's considered as column, the output is tt_tensor.
%   A and Y can be a {d,R} cell array. However, X is always a "single" TT
%   (no tensor chain), since that's how ALS works. Generally, X has the
%   same form as Y, except that it's always {d,1} in essense. X and Y can't
%   be sparse (SVD will destroy it anyways), but A can be.
%
%   Options are provided in form
%   'PropertyName1',PropertyValue1,'PropertyName2',PropertyValue2 and so
%   on. The parameters are set to default (in brackets in the following)
%   The list of option names and default values are:
%       o x0 - initial approximation to AY [rand with ranks of Y(:,1)]
%       o nswp - maximal number of sweeps [20]
%       o verb - verbosity level, 0-silent, 1-sweep info, 2-block info [1]
%       o kickrank - compression rank of the error, 
%                    i.e. enrichment size [4]
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
verb = 1;
init_qr = true;
X = [];

for i=1:2:length(varargin)-1
    switch lower(varargin{i})
        case 'nswp'
            nswp=varargin{i+1};
        case 'x0'
            X=varargin{i+1};
        case 'verb'
            verb=varargin{i+1};
        case 'kickrank'
            kickrank=varargin{i+1};
        case 'init_qr'
            init_qr = varargin{i+1};            
            
        otherwise
            warning('Unrecognized option: %s\n',varargin{i});
    end
end

% Grumble inputs
% First, the vector
if (isa(Y,'tt_tensor'))
    d = Y.d;
    m = Y.n;
    k = ones(d,1);
    ry = Y.r;
    Ry = 1;
    Y = core2cell(Y);
    outtype = 0;
elseif (isa(Y,'tt_matrix'))
    d = Y.d;
    m = Y.n;
    k = Y.m;
    ry = Y.r;
    Ry = 1;
    Y = core2cell(Y);
    outtype = 1;
else % {d,R}
    [d,Ry] = size(Y);
    m = zeros(d,1);
    k = zeros(d,1);
    ry = ones(d+1,Ry);
    for j=1:Ry
        for i=1:d
            [r1,n1,n2,r2]=size(Y{i,j});
            if (r1~=ry(i,j))
                error('Inconsistent ry(%d)', i);
            end;
            m(i) = n1;
            k(i) = n2;
            ry(i+1,j) = r2;
        end;
    end;
    outtype = 2;
end;

% Grumble matrix
if (isa(A, 'tt_matrix'))
    % Copy TT ranks and blocks from a tt_matrix
    ra = A.r;
    if (A.d~=d)||(any(A.m~=m))
        error('Mismatching TT dimensionalities in A and Y');
    end;
    n = A.n;
    Ra = 1;
    A = core2cell(A);
    for i=1:d
        A{i} = reshape(A{i}, ra(i)*n(i), m(i)*ra(i+1));
    end;
else
    % The matrix is given as a {d,R} cell array, meaning a sum
    % of R TT formats (needed for sparse blocks, in which case
    % it's the sum of R Kronecker products)
    Ra = size(A,2);
    if (size(A,1)~=d)
        error('Mismatching TT dimensionalities in A and Y');
    end;
    ra = ones(d+1,Ra);
    n = ones(d,1);
    for j=1:Ra
        for i=1:d
            if (issparse(A{i,j}))
                % Sparse matrix can only have TT ranks 1
                [n1,n2]=size(A{i,j});
                r2 = n2/m(i);
                n1 = n1/ra(i,j);
                if(abs(r2-round(r2))>sqrt(eps))
                    error('A is sparse, but the column size are not divisible by m');
                end;
                n2 = m(i);
            else
                [r1,n1,n2,r2]=size(A{i,j});
                if (r1~=ra(i,j))
                    error('Inconsistent ra(%d)', i);
                end;
                A{i,j} = reshape(A{i,j}, r1*n1, n2*r2);
            end;
            if (n2~=m(i))
                error('Mismatching m in A');
            end;
            n(i) = n1;
            ra(i+1,j) = r2;
        end;
    end;
end;

% Initial guess
if (isempty(X))
    X = tt_rand(n.*k, d, ry(:,1), 1);
    rx = X.r;
    X = core2cell(X);
    init_qr = false;
else
    % X was given, parse it
    if (isa(X,'tt_tensor'))
        if (X.d~=d)||(any(X.n~=n.*k))
            error('Mismatching dimensions in x0');
        end;
        rx = X.r;
        X = core2cell(X);
    elseif (isa(Y,'tt_matrix'))
        if (X.d~=d)||(any(X.n~=n))||(any(X.m~=k))
            error('Mismatching dimensions in x0');
        end;        
        rx = X.r;
        X = core2cell(X);
    else % {d,R}
        X = X(:,end); % There could be several TC columns
        rx = ones(d+1,1);
        for i=1:d
            [r1,n1,n2,r2]=size(X{i});
            if (r1~=rx(i))
                error('Inconsistent rx(%d)', i);
            end;
            if (n1~=n(i))||(n2~=k(i))
                error('Inconsistent n/k in x0');
            end;
            rx(i+1,j) = r2;
        end;
    end;
end;

% Reductions
XAY = cell(d+1,Ra,Ry);
ZAY = cell(d+1,Ra,Ry);
ZX = cell(d+1,1);
ZX{1} = 1; ZX{d+1} = 1;
for j=1:Ry
    for i=1:Ra
        XAY{1,i,j} = 1; XAY{d+1,i,j} = 1;
        ZAY{1,i,j} = 1; ZAY{d+1,i,j} = 1;
    end;
end;
% Residual rank
rz = [1; kickrank*ones(d-1,1); 1];

nrms = ones(d,1); % To prevent doulbe overflow

% Initial ort
for i=1:d-1
    if (init_qr)
        cr = reshape(X{i}, rx(i)*n(i)*k(i), rx(i+1));
        [cr,R]=qr(cr, 0);
        nrmr = norm(R, 'fro');
        if (nrmr>0)
            R = R/nrmr;
        end;        
        cr2 = reshape(X{i+1}, rx(i+1), n(i+1)*k(i+1)*rx(i+2));
        cr2 = R*cr2;
        rx(i+1) = size(cr, 2);
        X{i} = reshape(cr, rx(i), n(i), k(i), rx(i+1));
        X{i+1} = reshape(cr2, rx(i+1), n(i+1), k(i+1), rx(i+2));
    end;
    % Reduce
    [XAY(i+1,:,:),nrms(i)] = leftreduce_matrix(XAY(i,:,:), X{i}, A(i,:), Y(i,:), rx(i),n(i),k(i),rx(i+1), Ra,ra(i,:),ra(i+1,:), Ry,ry(i,:),m(i),ry(i+1,:));
    % Residual reductions    
    if (kickrank>0)
        cr = randn(rz(i)*n(i)*k(i), rz(i+1));
        [cr,~]=qr(cr, 0);
        rz(i+1) = size(cr, 2);        
        ZAY(i+1,:,:) = leftreduce_matrix(ZAY(i,:,:), cr, A(i,:), Y(i,:), rz(i),n(i),k(i),rz(i+1), Ra,ra(i,:),ra(i+1,:), Ry,ry(i,:),m(i),ry(i+1,:), nrms(i));
        ZX(i+1) = leftreduce_vector(ZX(i), cr, X(i), rz(i),n(i),k(i),rz(i+1), 1,rx(i),rx(i+1));
    end;    
end;

i = d;
dir = -1;
swp = 1;
max_dx = 0;
% Iteration
while (swp<=nswp)
    % Project the MatVec generating vector
    cr = local_matvec(Y(i,:), Ry,ry(i,:),m(i),k(i),ry(i+1,:), rx(i),n(i),rx(i+1), XAY(i,:,:), A(i,:), XAY(i+1,:,:), Ra,ra(i,:),ra(i+1,:));
    nrms(i) = norm(cr, 'fro');
    % The main goal is to keep y{i} of norm 1
    if (nrms(i)>0)
        cr = cr/nrms(i);
    else
        nrms(i)=1;
    end;    
    X{i} = reshape(X{i}, rx(i)*n(i)*k(i)*rx(i+1), 1);
    dx = norm(cr-X{i});
    max_dx = max(max_dx, dx);
    
    % Truncation and enrichment
    if ((dir>0)&&(i<d))
        cr = reshape(cr, rx(i)*n(i)*k(i), rx(i+1));
        [u,s,v]=svd(cr, 'econ');
        s = diag(s);
        r = my_chop2(s, tol*norm(s)/sqrt(d));
        u = u(:,1:r);
        v = conj(v(:,1:r))*diag(s(1:r));
        
        % Prepare enrichment, if needed
        if (kickrank>0)
            cr = u*v.';
            cr = reshape(cr, rx(i)*n(i)*k(i), rx(i+1));
            % For updating z            
            crz = local_matvec(Y(i,:), Ry,ry(i,:),m(i),k(i),ry(i+1,:), rz(i),n(i),rz(i+1), ZAY(i,:,:), A(i,:), ZAY(i+1,:,:), Ra,ra(i,:),ra(i+1,:));
            crz = reshape(crz, rz(i)*n(i)*k(i), rz(i+1));            
            ys = cr*ZX{i+1};
            yz = reshape(ys, rx(i), n(i)*k(i)*rz(i+1));
            yz = ZX{i}*yz;
            yz = reshape(yz, rz(i)*n(i)*k(i), rz(i+1));
            crz = crz/nrms(i) - yz;            
            nrmz = norm(crz,'fro');
            % For adding into solution
            crs = local_matvec(Y(i,:), Ry,ry(i,:),m(i),k(i),ry(i+1,:), rx(i),n(i),rz(i+1), XAY(i,:,:), A(i,:), ZAY(i+1,:,:), Ra,ra(i,:),ra(i+1,:));
            crs = reshape(crs, rx(i)*n(i)*k(i), rz(i+1));
            crs = crs/nrms(i) - ys;
            u = [u,crs];
            [u,R]=qr(u, 0);
            v = [v, zeros(rx(i+1), rz(i+1))];
            v = v*R.';
            r = size(u, 2);
        end;
        X{i} = reshape(u, rx(i), n(i), k(i), r);        
        cr2 = reshape(X{i+1}, rx(i+1), n(i+1)*k(i+1)*rx(i+2));
        v = reshape(v, rx(i+1), r);
        cr2 = v.'*cr2;
        X{i+1} = reshape(cr2, r, n(i+1), k(i+1), rx(i+2));
        
        rx(i+1) = r;
        
        nrms(i+dir) = nrms(i); % this is the norm of my block, save it
        % Reduce
        [XAY(i+1,:,:),nrms(i)] = leftreduce_matrix(XAY(i,:,:), X{i}, A(i,:), Y(i,:), rx(i),n(i),k(i),rx(i+1), Ra,ra(i,:),ra(i+1,:), Ry,ry(i,:),m(i),ry(i+1,:));
        % Enrichment
        if (kickrank>0)
            [crz,~]=qr(crz, 0);
            rz(i+1) = size(crz, 2);
            ZAY(i+1,:,:) = leftreduce_matrix(ZAY(i,:,:), crz, A(i,:), Y(i,:), rz(i),n(i),k(i),rz(i+1), Ra,ra(i,:),ra(i+1,:), Ry,ry(i,:),m(i),ry(i+1,:), nrms(i));
            ZX(i+1) = leftreduce_vector(ZX(i), crz, X(i), rz(i),n(i),k(i),rz(i+1), 1,rx(i),rx(i+1));
        end;
        
    elseif ((dir<0)&&(i>1))        
        cr = reshape(cr, rx(i), n(i)*k(i)*rx(i+1));
        [u,s,v]=svd(cr, 'econ');
        s = diag(s);
        r = my_chop2(s, tol*norm(s)/sqrt(d));
        v = conj(v(:,1:r));
        u = u(:,1:r)*diag(s(1:r));
        
        % Prepare enrichment, if needed
        if (kickrank>0)
            cr = u*v.';
            cr = reshape(cr, rx(i), n(i)*k(i)*rx(i+1));
            % For updating z
            crz = local_matvec(Y(i,:), Ry,ry(i,:),m(i),k(i),ry(i+1,:), rz(i),n(i),rz(i+1), ZAY(i,:,:), A(i,:), ZAY(i+1,:,:), Ra,ra(i,:),ra(i+1,:));
            crz = reshape(crz, rz(i), n(i)*k(i)*rz(i+1));            
            ys = ZX{i}*cr;
            yz = reshape(ys, rz(i)*n(i)*k(i), rx(i+1));
            yz = yz*ZX{i+1};
            yz = reshape(yz, rz(i), n(i)*k(i)*rz(i+1));
            crz = crz/nrms(i) - yz;
            nrmz = norm(crz,'fro');
            % For adding into solution
            crs = local_matvec(Y(i,:), Ry,ry(i,:),m(i),k(i),ry(i+1,:), rz(i),n(i),rx(i+1), ZAY(i,:,:), A(i,:), XAY(i+1,:,:), Ra,ra(i,:),ra(i+1,:));
            crs = reshape(crs, rz(i), n(i)*k(i)*rx(i+1));
            crs = crs/nrms(i) - ys;
            v = [v,crs.'];
            [v,R]=qr(v, 0);
            u = [u, zeros(rx(i), rz(i))];
            u = u*R.';
            r = size(v, 2);
        end;
        cr2 = reshape(X{i-1}, rx(i-1)*n(i-1)*k(i-1), rx(i));
        cr2 = cr2*u;
        X{i-1} = reshape(cr2, rx(i-1), n(i-1), k(i-1), r);
        X{i} = reshape(v.', r, n(i), k(i), rx(i+1));
        
        rx(i) = r;
        
        nrms(i+dir) = nrms(i); % this is the norm of my block, save it
        % Reduce
        [XAY(i,:,:),nrms(i)] = rightreduce_matrix(XAY(i+1,:,:), X{i}, A(i,:), Y(i,:), rx(i),n(i),k(i),rx(i+1), Ra,ra(i,:),ra(i+1,:), Ry,ry(i,:),m(i),ry(i+1,:));
        % Enrich
        if (kickrank>0)
            [crz,~]=qr(crz.', 0);
            rz(i) = size(crz, 2);            
            ZAY(i,:,:) = rightreduce_matrix(ZAY(i+1,:,:), crz, A(i,:), Y(i,:), rz(i),n(i),k(i),rz(i+1), Ra,ra(i,:),ra(i+1,:), Ry,ry(i,:),m(i),ry(i+1,:), nrms(i));
            ZX(i) = rightreduce_vector(ZX(i+1), crz, X(i), rz(i),n(i),k(i),rz(i+1), 1,rx(i),rx(i+1));
        end;
    else
        X{i} = reshape(cr, rx(i), n(i), k(i), rx(i+1));
    end;
    
    if (verb>1)
        fprintf('amen-mm: swp=[%d,%d], dx=%3.3e, r=%d, |X|=%3.3e, |z|=%3.3e\n', swp, i, dx, r, norm(cr,'fro'), nrmz);
    end;
    
    % Stopping or reversing
    if ((dir>0)&&(i==d))||((dir<0)&&(i==1))
        if (verb>0)
            fprintf('amen-mm: swp=%d{%d}, max_dx=%3.3e, max_r=%d\n', swp, (1-dir)/2, max_dx, max(rx));
        end;
        if ((max_dx<tol)||(swp==nswp))&&(dir>0)
            break;
        end;
        swp = swp+1;
        max_dx = 0;
        dir = -dir;
    else
        i = i+dir;
    end;
end;

% Distribute norms equally...
nrms = exp(sum(log(nrms))/d);
% ... and plug them into y
for i=1:d
    X{i} = X{i}*nrms;
end;

% Return the correct form
if (outtype==0)
    for i=1:d
        X{i} = reshape(X{i}, rx(i), n(i)*k(i), rx(i+1));        
    end;
    X = cell2core(tt_tensor, X);
elseif (outtype==1)
    X = cell2core(tt_matrix, X);
end;

end





% Accumulates the left reduction W{1:k}'*A{1:k}*X{1:k}
function [WAX2,nrm] = leftreduce_matrix(WAX2, w, A, x, rw1,n,k,rw2, Ra,ra1,ra2, Rx,rx1,m,rx2, extnrm)
% Left WAX has the form of the first matrix TT block, i.e. [rw, rx, ra]
if (nargin<16)
    extnrm = [];
end;
if (nargout>1)
    nrm = 0;
end;

w = reshape(w, rw1*n, k, rw2);
w = permute(w, [1,3,2]);
w = reshape(w, rw1, n*rw2*k);
for j=1:Rx
    xc = reshape(x{j}, rx1(j)*m, k, rx2(j));
    xc = permute(xc, [3,2,1]);
    xc = reshape(xc, rx2(j), k*rx1(j)*m);
    for i=1:Ra
        WAX2{1,i,j} = reshape(WAX2{1,i,j}, rw1, rx1(j)*ra1(i));
        WAX2{1,i,j} = w'*WAX2{1,i,j}; % size n rw2 x rx1 ra1
        WAX2{1,i,j} = reshape(WAX2{1,i,j}, n, rw2*k*rx1(j)*ra1(i));
        WAX2{1,i,j} = WAX2{1,i,j}.';
        WAX2{1,i,j} = reshape(WAX2{1,i,j}, rw2*k*rx1(j), ra1(i)*n);
        WAX2{1,i,j} = WAX2{1,i,j}*A{i}; % size rw2 rx1 m ra2
        WAX2{1,i,j} = reshape(WAX2{1,i,j}, rw2, k*rx1(j)*m*ra2(i));
        WAX2{1,i,j} = WAX2{1,i,j}.';
        WAX2{1,i,j} = reshape(WAX2{1,i,j}, k*rx1(j)*m, ra2(i)*rw2);
        WAX2{1,i,j} = xc*WAX2{1,i,j}; % size rx2, ra2 rw2
        WAX2{1,i,j} = reshape(WAX2{1,i,j}, rx2(j)*ra2(i), rw2);
        WAX2{1,i,j} = WAX2{1,i,j}.';
%         WAX2{i,j} = reshape(WAX2{i,j}, rw2, rx2(j), ra2(i));
        if (nargout>1)
            nrm = max(nrm, norm(WAX2{1,i,j}, 'fro'));
        end;
    end;
end;
if (nargout>1)
    % Extract the scale to prevent overload
    if (nrm>0)
        for i=1:Rx*Ra
            WAX2{i} = WAX2{i}/nrm;
        end;
    else
        nrm=1;
    end;
elseif (~isempty(extnrm))
    % Override the normalization
    for i=1:Rx*Ra
        WAX2{i} = WAX2{i}/extnrm;
    end;
end;
end

% Accumulates the left reduction W{1:k}'*X{1:k}
function [WX2] = leftreduce_vector(WX1, w, x, rw1,n,k,rw2, Rx,rx1,rx2, extnrm)
% Left WX has the form of the first vector TT block, i.e. [rw, rx]
if (nargin<11)
    extnrm = [];
end;

WX2 = WX1;
wc = reshape(w, rw1, n*k*rw2);
for i=1:Rx
    WX2{1,i} = wc'*WX2{1,i}; % size n rw2 x rx1
    WX2{1,i} = reshape(WX2{1,i}, n*k, rw2*rx1(i));
    WX2{1,i} = WX2{1,i}.';
    WX2{1,i} = reshape(WX2{1,i}, rw2, rx1(i)*n*k);
    tmp = reshape(x{i}, rx1(i)*n*k, rx2(i));
    WX2{1,i} = WX2{1,i}*tmp; % size rw2, rx2
end;
if (~isempty(extnrm))
    % Override the normalization
    for i=1:Rx
        WX2{i} = WX2{i}/extnrm;
    end;
end;
end

% Accumulates the right reduction W{k:d}'*A{k:d}*X{k:d}
function [WAX1,nrm] = rightreduce_matrix(WAX1, w, A, x, rw1,n,k,rw2, Ra,ra1,ra2, Rx,rx1,m,rx2, extnrm)
% Right WAX has the form of the last matrix TT block, i.e. [ra, rw, rx]
if (nargin<16)
    extnrm = [];
end;
if (nargout>1)
    nrm = 0;
end;

w = reshape(w, rw1*n, k, rw2);
w = permute(w, [1,3,2]);
w = reshape(w, rw1, n*rw2*k);
w = conj(w);
for j=1:Rx
    xc = reshape(x{j}, rx1(j)*m, k, rx2(j));
    xc = permute(xc, [2,1,3]);
    xc = reshape(xc, k*rx1(j)*m, rx2(j));
    for i=1:Ra
        WAX1{1,i,j} = reshape(WAX1{1,i,j}, ra2(i)*rw2, rx2(j));
        WAX1{1,i,j} = xc*WAX1{1,i,j}.'; % size rx1 m x ra2 rw2
        WAX1{1,i,j} = reshape(WAX1{1,i,j}, k*rx1(j), m*ra2(i)*rw2);
        WAX1{1,i,j} = WAX1{1,i,j}.';
        WAX1{1,i,j} = reshape(WAX1{1,i,j}, m*ra2(i), rw2*k*rx1(j));
        WAX1{1,i,j} = A{i}*WAX1{1,i,j}; % size ra1(k)*n, rw2*rx1
        WAX1{1,i,j} = reshape(WAX1{1,i,j}, ra1(i), n*rw2*k*rx1(j));
        WAX1{1,i,j} = WAX1{1,i,j}.';
        WAX1{1,i,j} = reshape(WAX1{1,i,j}, n*rw2*k, rx1(j)*ra1(i));
        WAX1{1,i,j} = w*WAX1{1,i,j}; % size rw1, rx1 ra1
        WAX1{1,i,j} = reshape(WAX1{1,i,j}, rw1*rx1(j), ra1(i));
        WAX1{1,i,j} = WAX1{1,i,j}.';
%         WAX1{i,j} = reshape(WAX1{i,j}, ra1(i), rw1, rx1(j));
        if (nargout>1)
            nrm = max(nrm, norm(WAX1{1,i,j}, 'fro'));
        end;
    end;
end;
if (nargout>1)
    % Extract the scale to prevent overload
    if (nrm>0)
        for i=1:Rx*Ra
            WAX1{i} = WAX1{i}/nrm;
        end;
    else
        nrm=1;
    end;
elseif (~isempty(extnrm))
    % Override the normalization
    for i=1:Rx*Ra
        WAX1{i} = WAX1{i}/extnrm;
    end;
end;
end


% Accumulates the right reduction W{k:d}'*X{k:d}
function [WX1] = rightreduce_vector(WX2, w, x, rw1,n,k,rw2, Rx,rx1,rx2, extnrm)
% Right WX has the form of the last vector TT block, i.e. [rx, rw]
if (nargin<11)
    extnrm = [];
end;

WX1 = WX2;
wc = reshape(w, rw1, n*k*rw2);
for i=1:Rx
    tmp = reshape(x{i}, rx1(i)*n*k, rx2(i));
    WX1{1,i} = tmp*WX1{1,i}; % size rx1 n x rw2
    WX1{1,i} = reshape(WX1{1,i}, rx1(i), n*k*rw2);
    WX1{1,i} = WX1{1,i}*wc'; % size rx1, rw1
end;
if (~isempty(extnrm))
    % Override the normalization
    for i=1:Rx
        WX1{i} = WX1{i}/extnrm;
    end;
end;
end


% A matrix-matrix product for the matrix in the 3D TT (WAX1-A-WAX2), and
% full matrix of size (rx1*m*k*rx2). Returns (rw1*n*k*rw2)
function [w]=local_matvec(x, Rx,rx1,m,k,rx2, rw1,n,rw2, WAX1, A, WAX2, Ra,ra1,ra2)
w = zeros(rw1*n*rw2,k);
for j=1:Rx
    xc = reshape(x{j}, rx1(j)*m, k, rx2(j));
    xc = permute(xc, [2,1,3]);
    xc = reshape(xc, k*rx1(j)*m, rx2(j));
    for i=1:Ra
        tmp = reshape(WAX2{1,i,j}, ra2(i)*rw2, rx2(j));
        wk = xc*tmp.';
        wk = reshape(wk, k*rx1(j), m*ra2(i)*rw2);
        wk = wk.';
        wk = reshape(wk, m*ra2(i), rw2*k*rx1(j));
        wk = A{i}*wk;
        wk = reshape(wk, ra1(i)*n*rw2*k, rx1(j));
        wk = wk.';
        wk = reshape(wk, rx1(j)*ra1(i), n*rw2*k);
        tmp = reshape(WAX1{1,i,j}, rw1, rx1(j)*ra1(i));
        wk = tmp*wk;
        wk = reshape(wk, rw1*n*rw2, k);
        w = w+wk;
    end;
end;
w = reshape(w, rw1*n, rw2, k);
w = permute(w, [1,3,2]);
w = reshape(w,[],1);
end
