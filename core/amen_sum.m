function [y,z]=amen_sum(X, c, tol, varargin)
%Approximate the linear combination of TT tensors via the AMEn iteration
%   [y,z]=amen_sum(X, v, tol, varargin)
%   Attempts to approximate the y(j) = sum_i X{i}*c(i,j)
%   with accuracy TOL using the AMEn+ALS iteration.
%
%   Options are provided in form
%   'PropertyName1',PropertyValue1,'PropertyName2',PropertyValue2 and so
%   on. The parameters are set to default (in brackets in the following)
%   The list of option names and default values are:
%       o y0 - initial approximation to y [rand rank-2]
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
%       o can - whether the input is in CP format [false]
%       o multrank - shall we try to grow ranks multiplicatively by adding
%       the previous iterand [false]
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

% Special case: amen_sum called from another amen package
if (~isempty(varargin))
    v1 = varargin{1};
    if (isa(v1, 'cell'))
        varargin=v1;
    end;
end;

% Default valuse
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
rmax = Inf;
multrank = false;

can = false; % Whether the input is the canonical format

% Read parameters from input sequence
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
        case 'rmax'
            rmax = varargin{i+1};
        case 'can'
            can = varargin{i+1};
        case 'multrank'
            multrank = varargin{i+1};            
            
        otherwise
            warning('Unrecognized option: %s\n',varargin{i});
    end
end

% Convert input to canonical or TT cell array
N = size(c,1);
M = size(c,2);
Xin = X;

if (can)
    N = 1;
    % X itself is already given in the required format
    d = numel(X);
    R = size(X{1}, 2);    
    nrms = ones(d,1);
    n = zeros(d,1);
    for i=1:d
        n(i) = size(X{i}, 1);
    end;
    vectype = 1;
else
    R = 1;
    for j=1:N
        if (isa(Xin{j}, 'tt_tensor'))
            if (j==1)
                d = Xin{j}.d;
                n = Xin{j}.n;
                X = cell(d,N);
                
                % A storage for norms.
                % Z may be kept normalized, but for y we will have y_real = y*prod(nrms).
                % The same for all x. But note that we keep only the _scale_ for all x vecs.
                nrms = ones(d,1);
            end;
            X(:,j) = core2cell(Xin{j});
            vectype = 1; % tt_tensor
        else
            if (j==1)
                d = numel(Xin{j});
                n = zeros(d,1);
                for i=1:d
                    n(i) = size(Xin{j}{i}, 2);
                end;
                X = cell(d,N);
            end;
            X(:,j) = Xin{j};
            vectype = 0; % cell
        end;
    end;
end;

% Initial guess
if (isempty(y))
    init_qr = false;
    [y,ry] = gen_rand(n,d,2);
else
    if (isa(y, 'tt_tensor'))
        y = core2cell(y);
    end;
    ry = ones(d+1,1);
    for i=1:d
        ry(i+1) = size(y{i}, 3);
    end;
end;

% Enrichment vector
if (kickrank+kickrank2>0)
    if (isempty(z))
        [z,rz] = gen_rand(n,d,kickrank+kickrank2);
        init_qr_z = false;
    else
        init_qr_z = true;
        if (isa(z, 'tt_tensor'))
            z = core2cell(z);
        end;
        rz = ones(d+1,1);
        for i=1:d
            rz(i+1) = size(z{i}, 3);
        end;
    end;
    
    phizx = cell(d+1,N);
    for j=1:N
        phizx{1,j}=ones(1,R); phizx{d+1,j}=ones(R,1);
    end;
    phizy = cell(d+1,1); phizy{1}=1; phizy{d+1}=1;
end;

phiyx = cell(d+1,N);
for j=1:N
    phiyx{1,j}=ones(1,R); phiyx{d+1,j}=ones(R,1);
end;

% Initial ort
for i=1:d-1
    if (init_qr)
        cr = reshape(y{i}, ry(i)*n(i), ry(i+1));
        if (strcmp(renorm, 'gram'))&&(ry(i)*n(i)>5*ry(i+1))
            [cr,~,R]=svdgram(cr);
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
    % We will extract the norms of Y'X
    [phiyx(i+1,:),nrms(i)] = compute_next_Phi(phiyx(i,:), y{i}, X(i,:), 'lr', can);
    
    if (kickrank+kickrank2>0)
        if (init_qr_z)
            cr = reshape(z{i}, rz(i)*n(i), rz(i+1));
            if (strcmp(renorm, 'gram'))&&(rz(i)*n(i)>5*rz(i+1))
                [cr,~,R]=svdgram(cr);
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
        end;
        
        % But we have to renorm Z'X/|Y'X| to keep Z in the same
        % scale. Hope it will not make |Z'X| too small
        phizx(i+1,:) = compute_next_Phi(phizx(i,:), z{i}, X(i,:), 'lr', can, nrms(i));
        % Z'Y is actually an orthogonal matrix and does not require renorm
        phizy(i+1) = compute_next_Phi(phizy(i), z{i}, y(i), 'lr');
    end;
end;

i = d;
dir = -1;
swp = 1;
max_dx = 0;
y{d} = reshape(y{d}, ry(d)*n(d), ry(d+1));
y{d} = y{d}(:,1:min(ry(d+1), M));
y{d} = [y{d}, zeros(ry(d)*n(d), M-ry(d+1))];
y{d} = reshape(y{d}, ry(d), n(d), 1, M);
ry(d+1) = 1;

while (swp<=nswp)
    % Project the sum
    cry = proj_sum(phiyx(i,:), X(i,:), phiyx(i+1,:), c);
    % It will survive in the terminal block. All we need.
    nrms(i) = norm(cry, 'fro');
    % The main goal is to keep y{i} of norm 1
    if (nrms(i)>0)
        cry = cry/nrms(i);
    else
        nrms(i)=1;
    end;
    
    % Check stopping criteria
    y{i} = reshape(y{i}, ry(i)*n(i)*ry(i+1), M);
    dx = norm(cry-y{i}, 'fro')/norm(cry, 'fro');
    max_dx = max(max_dx, dx);
    
    % Truncation and enrichment
    if ((dir>0)&&(i<d))
        cry = reshape(cry, ry(i)*n(i), ry(i+1)*M);
        if (M==1)&&(multrank)
            % Try to accelerate the rank growth
            cry = [cry, reshape(y{i}, ry(i)*n(i), ry(i+1))];
        end;
        if (strcmp(renorm, 'gram'))
            [u,~,v]=svdgram(cry, tol/sqrt(d));
            v = v.';
            if (size(u,2)>rmax)
                u = u(:,1:rmax);
                v = v(:,1:rmax);
            end;
            r = size(u,2);
        else
            [u,s,v]=svd(cry, 'econ');
            s = diag(s);
            r = my_chop2(s, tol*norm(s)/sqrt(d));
            r = min(r,rmax);
            u = u(:,1:r);
            v = conj(v(:,1:r))*diag(s(1:r));
        end;        
        if (M==1)&&(multrank)
            % Remove auxiliary previous solution
            v = reshape(v, ry(i+1), 2, r);
            v = v(:,1,:);
            v = reshape(v, ry(i+1), r);
        end;
        
        % Prepare enrichment, if needed
        if (kickrank+kickrank2>0)
            cry = u*v.';
            cry = reshape(cry, ry(i)*n(i)*ry(i+1), M);
            % For updating z
            crys = cry.';
            crys = reshape(crys, M*ry(i)*n(i), ry(i+1));
            crys = crys*phizy{i+1};
            crys = reshape(crys, M, ry(i)*n(i)*rz(i+1));
            crys = crys.';
            cryz = reshape(crys, ry(i), n(i)*rz(i+1)*M);
            cryz = phizy{i}*cryz;
            cryz = reshape(cryz, rz(i)*n(i)*rz(i+1), M);
            
            crz = proj_sum(phizx(i,:), X(i,:), phizx(i+1,:), c);
            crz = crz/nrms(i) - cryz;
            nrmz = norm(crz, 'fro');
            crz = reshape(crz, rz(i)*n(i), rz(i+1)*M);
            [crz,~,~]=svd(crz, 'econ');
            crz = crz(:, 1:min(size(crz,2), kickrank));
            if (kickrank2>0)
                crz = [crz, randn(rz(i)*n(i), kickrank2)];
                [crz,~]=qr(crz,0);
            end;
            % For adding into solution
            if (fkick)
                crs = proj_sum(phiyx(i,:), X(i,:), phizx(i+1,:), c);
                crs = crs/nrms(i) - crys;
                crs = reshape(crs, ry(i)*n(i), rz(i+1)*M);
                [crs,~,~]=svd(crs, 'econ');
                crs = crs(:, 1:min(size(crs,2), kickrank));
                u = [u,crs];
                if (strcmp(renorm, 'gram'))&&(ry(i)*n(i)>5*(ry(i+1)+rz(i+1)))
                    [u,~,R]=svdgram(u);
                else
                    [u,R]=qr(u, 0);
                end;
                v = [v, zeros(ry(i+1)*M, size(crs,2))];
                v = v*R.';
                r = size(u, 2);
            end;
        end;
        y{i} = reshape(u, ry(i), n(i), r);
        
        cr2 = reshape(y{i+1}, ry(i+1), n(i+1)*ry(i+2));
        v = reshape(v, ry(i+1), M*r);
        cr2 = v.'*cr2;
        cr2 = reshape(cr2, M, r*n(i+1)*ry(i+2));
        y{i+1} = reshape(cr2.', r, n(i+1), ry(i+2), M);
        
        ry(i+1) = r;
        
        [phiyx(i+1,:), nrms(i)] = compute_next_Phi(phiyx(i,:), y{i}, X(i,:), 'lr', can);
        
        if (kickrank+kickrank2>0)
            rz(i+1) = size(crz, 2);
            z{i} = reshape(crz, rz(i), n(i), rz(i+1));
            % z{i+1} will be recomputed from scratch in the next step
            phizx(i+1,:) = compute_next_Phi(phizx(i,:), z{i}, X(i,:), 'lr', can, nrms(i));
            phizy(i+1) = compute_next_Phi(phizy(i), z{i}, y(i), 'lr');
        end;
    elseif ((dir<0)&&(i>1))
        cry = reshape(cry.', M*ry(i), n(i)*ry(i+1));
        if (M==1)&&(multrank)
            % Try to accelerate the rank growth
            cry = [cry; reshape(y{i}, ry(i), n(i)*ry(i+1))];
        end;
        [u,s,v]=svd(cry, 'econ');
        s = diag(s);
        r = my_chop2(s, tol*norm(s)/sqrt(d));
        r = min(r,rmax);
        u = u(:,1:r)*diag(s(1:r));
        v = conj(v(:,1:r));
        
        if (M==1)&&(multrank)
            % Remove auxiliary previous solution
            u = reshape(u, ry(i), 2, r);
            u = u(:,1,:);
            u = reshape(u, ry(i), r);
        end;             
        
        % Prepare enrichment, if needed
        if (kickrank+kickrank2>0)
            cry = u*v.';
            cry = reshape(cry, M, ry(i)*n(i)*ry(i+1));
            % For updating z
            crys = cry.';
            crys = reshape(crys, ry(i), n(i)*ry(i+1)*M);
            crys = phizy{i}*crys;
            crys = reshape(crys, rz(i)*n(i)*ry(i+1), M);
            cryz = reshape(crys.', M*rz(i)*n(i), ry(i+1));
            cryz = cryz*phizy{i+1};
            cryz = reshape(cryz, M, rz(i)*n(i)*rz(i+1));
            cryz = cryz.';
            crz = proj_sum(phizx(i,:), X(i,:), phizx(i+1,:), c);
            crz = crz/nrms(i) - cryz;
            nrmz = norm(crz, 'fro');
            crz = reshape(crz.', M*rz(i), n(i)*rz(i+1));
            [~,~,crz]=svd(crz, 'econ');
            crz = crz(:, 1:min(size(crz,2), kickrank));
            if (kickrank2>0)
                crz = [crz, randn(n(i)*rz(i+1), kickrank2)];
                [crz,~]=qr(crz,0);
            end;            
            crz = crz';
            % To add into solution
            crs = proj_sum(phizx(i,:), X(i,:), phiyx(i+1,:), c);
            crs = crs/nrms(i) - crys;
            crs = reshape(crs.', M*rz(i), n(i)*ry(i+1));
            [~,~,crs]=svd(crs, 'econ');
            crs = crs(:, 1:min(size(crs,2), kickrank));
            v = [v,conj(crs)];
            if (strcmp(renorm, 'gram'))&&(ry(i+1)*n(i)>5*(ry(i)+rz(i)))
                [v,~,R]=svdgram(v);
            else
                [v,R]=qr(v, 0);
            end;
            u = [u, zeros(M*ry(i), size(crs,2))];
            u = u*R.';
            r = size(v, 2);
        end;
        y{i} = reshape(v.', r, n(i), ry(i+1));
        
        cr2 = reshape(y{i-1}, ry(i-1)*n(i-1), ry(i));
        u = reshape(u, M, ry(i)*r);
        u = u.';
        u = reshape(u, ry(i), r*M);
        cr2 = cr2*u;
        cr2 = reshape(cr2, ry(i-1), n(i-1), r, M);
        y{i-1} = cr2;
        
        ry(i) = r;
        
        [phiyx(i,:), nrms(i)] = compute_next_Phi(phiyx(i+1,:), y{i}, X(i,:), 'rl', can);
        
        if (kickrank+kickrank2>0)
            rz(i) = size(crz, 1);
            z{i} = reshape(crz, rz(i), n(i), rz(i+1));
            % z{i+1} will be recomputed from scratch in the next step
            phizx(i,:) = compute_next_Phi(phizx(i+1,:), z{i}, X(i,:), 'rl', can, nrms(i));
            phizy(i) = compute_next_Phi(phizy(i+1), z{i}, y(i), 'rl');
        end;
    else
        y{i} = reshape(y{i}, ry(i), n(i), ry(i+1), M);
    end;
    
    if (verb>1)
        fprintf('amen-sum: swp=[%d,%d], dx=%3.3e, r=%d, |y|=%3.3e, |z|=%3.3e\n', swp, i, dx, r, norm(cry, 'fro'), nrmz);
    end;
    
    % Stopping or reversing
    if ((dir>0)&&(i==d))||((dir<0)&&(i==1))
        if (verb>0)
            fprintf('amen-sum: swp=%d{%d}, max_dx=%3.3e, max_r=%d\n', swp, (1-dir)/2, max_dx, max(ry));
        end;
        if ((max_dx<tol)||(swp==nswp))&&(dir>0)
            break;
        end;
        max_dx = 0;
        if (dir>0); swp = swp+1; end;
        dir = -dir;
    else
        i = i+dir;
    end;
end;

y{d} = reshape(cry, ry(d), n(d), M);

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


function [Phi,nrm] = compute_next_Phi(Phi_prev, x, Y, direction, can, extnrm)
% Performs the recurrent Phi (or Psi) matrix computation
% Phi = Phi_prev * (x'y).
% If direction is 'lr', computes Psi
% if direction is 'rl', computes Phi
% A can be empty, then only x'y is computed.

% Phi1:  rx1, ry1
% Phi2:  ry2, rx2

if (nargin<5)||(isempty(can))
    can = false;
end;
if (nargin<6)
    extnrm = [];
end;

rx1 = size(x,1); n = size(x,2); rx2 = size(x,3);
N = numel(Y);
if (can) % We are working with the canonical format
    Phi = cell(1,1);
    R = size(Y{1},2);
    % Phi1: rx1, R
    % Phi2: R, rx2
    
    if (strcmp(direction, 'lr'))
        %lr: Phi1
        x = reshape(x, rx1, n*rx2);
        Phi{1} = x'*Phi_prev{1};
        Phi{1} = reshape(Phi{1}, n*rx2, R);
        Y = repmat(Y{1}, rx2, 1); % equalize the sizes for Hadamard product
        Phi{1} = Phi{1}.*Y;
        Phi{1} = reshape(Phi{1}, n, rx2*R);
        Phi{1} = sum(Phi{1}, 1);
        Phi{1} = reshape(Phi{1}, rx2, R);
        if (nargout>1)
            % Extract the scale to prevent overload
            nrm = norm(Phi{1}, 'fro');
            if (nrm>0)
                Phi{1} = Phi{1}/nrm;
            else
                nrm=1;
            end;
        elseif (~isempty(extnrm))
            % Override the normalization by the external one
            Phi{1} = Phi{1}/extnrm;
        end;
    else
        %rl: Phi2
        x = reshape(x, rx1*n, rx2);
        Phi{1} = Phi_prev{1}*x';
        Phi{1} = reshape(Phi{1}, R, rx1*n);
        Y = Y{1}.';
        Y = repmat(Y,rx1,1); % equalize the sizes for Hadamard product
        Y = reshape(Y, R, rx1*n);
        Phi{1} = Phi{1}.*Y;
        Phi{1} = reshape(Phi{1}, R*rx1, n);
        Phi{1} = sum(Phi{1}, 2);
        Phi{1} = reshape(Phi{1}, R, rx1);
        if (nargout>1)
            % Extract the scale to prevent overload
            nrm = norm(Phi{1}, 'fro');
            if (nrm>0)
                Phi{1} = Phi{1}/nrm;
            else
                nrm=1;
            end;
        elseif (~isempty(extnrm))
            % Override the normalization by the external one
            Phi{1} = Phi{1}/extnrm;
        end;
    end;
else % a set of TT-tensors
    Phi = cell(1,N);
    if (strcmp(direction, 'lr'))
        %lr: Phi1
        x = reshape(x, rx1, n*rx2);
        if (nargout>1)
            nrm = 0;
        end;
        for i=1:N
            y = Y{i};
            ry1 = size(y,1); ry2 = size(y,3);
            Phi{i} = x'*Phi_prev{i};
            Phi{i} = reshape(Phi{i}, n, rx2*ry1);
            Phi{i} = Phi{i}.';
            Phi{i} = reshape(Phi{i}, rx2, ry1*n);
            y = reshape(y, ry1*n, ry2);
            Phi{i} = Phi{i}*y;
            Phi{i} = reshape(Phi{i}, rx2, ry2);
            if (nargout>1)
                nrm = max(nrm, norm(Phi{i}, 'fro'));
            end;
        end;
        if (nargout>1)
            % Extract the scale to prevent overload
            if (nrm>0)
                for i=1:N
                    Phi{i} = Phi{i}/nrm;
                end;
            else
                nrm=1;
            end;
        elseif (~isempty(extnrm))
            % Override the normalization
            for i=1:N
                Phi{i} = Phi{i}/extnrm;
            end;
        end;
    else
        %rl: Phi2
        x = reshape(x, rx1, n*rx2);
        if (nargout>1)
            nrm = 0;
        end;
        for i=1:N
            y = Y{i};
            ry1 = size(y,1); ry2 = size(y,3);
            y = reshape(y, ry1*n, ry2);
            Phi{i} = y*Phi_prev{i};
            Phi{i} = reshape(Phi{i}, ry1, n*rx2);
            Phi{i} = Phi{i}*x';
            Phi{i} = reshape(Phi{i}, ry1, rx1);
            if (nargout>1)
                nrm = max(nrm, norm(Phi{i}, 'fro'));
            end;
        end;
        if (nargout>1)
            % Extract the scale to prevent overload
            if (nrm>0)
                for i=1:N
                    Phi{i} = Phi{i}/nrm;
                end;
            else
                nrm=1;
            end;
        elseif (~isempty(extnrm))
            % Override the normalization
            for i=1:N
                Phi{i} = Phi{i}/extnrm;
            end;
        end;
    end;
end;

end


function [x,r]=gen_rand(n,d,r)
% Generate an orthogonal random vector
if (numel(r)==1)
    r = [1; r*ones(d-1,1); 1];
end;
x = cell(d,1);
for i=1:d
    cr = randn(r(i)*n(i), r(i+1));
    [cr,~]=qr(cr,0);
    r(i+1) = size(cr,2);
    x{i} = reshape(cr, r(i), n(i), r(i+1));
end;
end

function [y]=proj_sum(Phi1, X, Phi2, c)
N = size(c,1);
M = size(c,2);

ry1 = size(Phi1{1},1);
ry2 = size(Phi2{1},2);

if (numel(X)==1) % Canonical format
    X = X{1};    
    n = size(X,1);
    Phi1 = repmat(Phi1{1}, n, 1); % size ry1*n, R
    X = reshape(X,  1, n*N);
    X = repmat(X, ry1, 1);
    X = reshape(X, ry1*n, N);
    y = X.*Phi1;
    c = repmat(c, ry2, 1); % size N*ry2, M
    c = reshape(c, N, ry2*M);
    Phi2 = repmat(Phi2{1}, 1, M); % size N, ry2*M
    Phi2 = Phi2.*c;
    y = y*Phi2;
    y = reshape(y, ry1*n*ry2, M);
    if (issparse(y))
        y = full(y);
    end;
else
    n = size(X{1}, 2);
    y = zeros(ry1*n*ry2, M);
    for j=1:N
        rx1 = size(X{j},1); rx2 = size(X{j},3);
        cry = reshape(X{j}, rx1, n*rx2);
        cry = Phi1{j}*cry;
        cry = reshape(cry, ry1*n, rx2);
        cry = cry*Phi2{j};
        cry = reshape(cry, ry1*n*ry2, 1);
        cry = cry*c(j,:);
        y = y+cry;
    end;
end;
end