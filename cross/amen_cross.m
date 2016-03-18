function [y,Y]=amen_cross(inp, fun, tol, varargin)
% Block cross with error-based enrichment.
% function [y]=amen_cross(inp, fun, tol, varargin)
% Tries to interpolate the function(s) via the error-enriched maxvol-cross.
% This version accepts both index-defined (traditional) functions, 
% and elementwise functions depending on other tt_tensors.
%
% If inp==n is a column vector of mode sizes:
%       fun = @(ind)fun(ind) is a sought function of index ind,
%       ind comes as an array of sizes M x d,
%       fun should return the block M x B (vectorized variant).
%       M=1 if vec=false, hence the return may be either 1 x B or B x 1.
% If inp=={x1,x2,...,xN} a cell array of tt_tensors x1,x2,...xN:
%       fun = (x)fun(x) is a sought function of elements of x=(x1,x2,...).
%       it should receive a 2d array V of sizes M x N, where the
%       first dimension stays for the reduced set of spatial indices, and the
%       second is the enumerator of vectors in X.
%       The returned sizes should be M x B, where B is the number of
%       components in FUNS.
% In addition to the obligatory first inputs, the second function of another type may be
% added via optional parameters 'auxinp', 'auxfun' (see below).
%
% Optional arguments are provided in the form
% 'PropertyName1',PropertyValue1,'PropertyName2',PropertyValue2 and so on. 
% The list of option names and default values:
%       o y0 - initial approximation [random rank-2 tensor]
%       o nswp - maximal number of sweeps [20]
%       o zrank - rank of the error approx. Z [5]
%       o zrank2 - size of the secondary random enrichment for Z [zrank/2]
%       o kickrank - actual enrichment size [2]
%       o tol_exit - stopping difference between consecutive iterations [tol]
%       o verb - verbosity level, 0-silent, 1-sweep info, 2-block info [1]
%       o vec - whether index-fun can accept and return vectorized values [false]
%       o exitdir - if 1, return after the forward sweep, if -1, return the
%                   backward sweep, if 0, after any [1]
%       o auxinp - secondary input data
%       o auxfun - secondary input function
%       o max_err_jumps - stop when error increased max_err_jumps times.
%       o trunc_method - truncate ranks using TT-SVD (svd) or (cross) over
%                        rn x r local samples. The latter seems more robust
%
%********
%   References for the alternating optimization with enrichment (AMEn):
%   S. Dolgov, D. Savostyanov.
%   http://arxiv.org/abs/1301.6068 
%   http://arxiv.org/abs/1304.1222
%   
%   Please send feedback to: {sergey.v.dolgov,dmitry.savostyanov}@gmail.com
%
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
vars = varargin;

y = [];
nswp = 20;
zrank = 5;
zrank2 = [];
kickrank = 2;
tol_exit = tol;
verb = 1;
vec = false;
exitdir=1;
max_err_jumps = 2;
trunc_method = 'cross';
% trunc_method = 'svd';

auxinp = [];
auxfun = [];

i = 1;
while (i<length(vars))
    switch lower(vars{i})
        case 'y0'
            y=vars{i+1};
        case 'nswp'
            nswp=vars{i+1};
        case 'kickrank'
            kickrank=vars{i+1};            
        case 'zrank'
            zrank=vars{i+1};
        case 'zrank2'
            zrank2=vars{i+1};            
        case 'tol_exit'
            tol_exit=vars{i+1};  
        case 'verb'
            verb = vars{i+1};      
        case 'vec'
            vec = vars{i+1};   
        case 'auxinp'
            auxinp = vars{i+1}; 
        case 'auxfun'
            auxfun = vars{i+1};             
        case 'exitdir'
            exitdir=vars{i+1};
        case 'max_err_jumps'
            max_err_jumps=vars{i+1};
        case 'trunc_method'
            trunc_method=vars{i+1};
        otherwise
            warning('Option %s was not recognized', vars{i});
    end;
    i=i+2;
end;

% Half of Z rank will be enriched randomly
if (isempty(zrank2))
    zrank2 = round(zrank/2);
end;

X = [];
ifun = [];
ffun = [];
n = [];
% Distinguish general and TT-fun inputs
if (isa(inp, 'cell'))
    % First input is an array of TT tensors
    X = inp;
    ffun = fun;
else
    ifun = fun;
    n = inp;
end;
if (~isempty(auxinp))&&(~isempty(auxfun))
    if (isa(auxinp, 'cell'))
        % Second input is an array of TT tensors
        if (isempty(ffun))
            X = auxinp;
            ffun = auxfun;
        else
            error('Cannot use ffun on both inputs');
        end;
    else
        if (isempty(ifun))
            ifun = auxfun;
            n = auxinp;
        else
            error('Cannot use ifun on both inputs');
        end;
    end;
end;

% If there is a TT-fun part, prepare it for computations
if (~isempty(X))
    nx = numel(X);
    d = X{1}.d;
    n = X{1}.n;
    rx = zeros(d+1,nx);
    crX = cell(d,nx);
    for i=1:nx
        rx(:,i) = X{i}.r;
        crX(:,i) = core2cell(X{i});
    end;
    % Interface matrices
    phiyx = cell(d+1,nx);
    % these are for stuffing z indices into x
    phizx = cell(d+1,nx);
    for i=1:nx
        phiyx{1,i}=1; phiyx{d+1,i}=1;
        phizx{1,i}=1; phizx{d+1,i}=1;
    end;
end;

d = numel(n);
if (isempty(y))
    y = tt_rand(n, d, 2, -1);
end;
ry = y.r;
y = core2cell(y);

if (zrank>0)
    rz = [1; zrank*ones(d-1,1); 1];
%     z = tt_rand(n,d,zrank,-1);
%     rz = z.r;
%     z = core2cell(z);
    % Interp matrices for z and y->z    
%     phizz = cell(d+1,1);
%     phizz{1} = 1;
%     phizz{d+1} = 1;
    phizy = cell(d+1,1);
    phizy{1} = 1;
    phizy{d+1} = 1;
end;

% Interp matrices for y
if (strcmp(trunc_method, 'svd'))
    phiyy = cell(d+1,1);
    phiyy{1} = 1;
    phiyy{d+1}=1;
end;

if (~isempty(ifun))
    % Some place to store global indices
    Jy = cell(d+1,1);
    Jz = cell(d+1,1);
end;

% Store factorized norm
nrms = ones(d,1);

% Initial orthog and indices: maxvol over y, z
for i=d:-1:2
    cry = reshape(y{i}, ry(i), n(i)*ry(i+1));
    [cry,rv]=qr(cry.', 0);
    rv = rv.';
    cry = cry.';
    if (~strcmp(trunc_method, 'svd'))
        ind = maxvol2(cry.');
        yy = cry(:,ind);
        cry = yy\cry;
        rv = rv*yy;
        nrms(i) = 1;
    end;
    cr2 = reshape(y{i-1}, ry(i-1)*n(i-1), ry(i));
    cr2 = cr2*rv;
    ry(i) = size(cry,1);
    y{i} = reshape(cry, ry(i), n(i), ry(i+1));
    y{i-1} = reshape(cr2, ry(i-1), n(i-1), ry(i));   
    
    if (strcmp(trunc_method, 'svd'))
        cry = reshape(y{i}, ry(i)*n(i), ry(i+1));
        cry = cry*phiyy{i+1};
        cry = reshape(cry, ry(i), n(i)*ry(i+1));
        ind = maxvol2(cry.');
        phiyy{i} = cry(:,ind);
        % Extract the scale
        nrms(i) = 1./min(svd(phiyy{i}));
        phiyy{i} = phiyy{i}.*nrms(i);
    end;
    if (~isempty(ifun))
        % Index function is present
        Jy{i} = indexmerge((1:n(i))', Jy{i+1});
        Jy{i} = Jy{i}(ind,:);
    end;
    if (~isempty(ffun))
        % Elem function is present
        % Interface matrices for TT-fun
        for j=1:nx
            phiyx{i,j} = reshape(crX{i,j}, rx(i,j)*n(i), rx(i+1,j));
            phiyx{i,j} = phiyx{i,j}*phiyx{i+1,j};
            phiyx{i,j} = reshape(phiyx{i,j}, rx(i,j), n(i)*ry(i+1));
            phiyx{i,j} = phiyx{i,j}(:, ind);
        end;
    end;
    
    if (zrank>0)
        % The same for residual
        crz = randn(n(i)*rz(i+1), rz(i));
%         crz = reshape(z{i}, rz(i), n(i)*rz(i+1));
        [crz,~]=qr(crz, 0);
        ind = maxvol2(crz);        
%         cr2 = reshape(z{i-1}, rz(i-1)*n(i-1), rz(i));
%         cr2 = cr2*rv.'*crz(ind,:).';
        rz(i) = size(crz,2);
%         z{i} = reshape(crz.', rz(i), n(i), rz(i+1));
%         z{i-1} = reshape(cr2, rz(i-1), n(i-1), rz(i));
        
%         crz = reshape(z{i}, rz(i)*n(i), rz(i+1));
%         crz = crz*phizz{i+1};
%         crz = reshape(crz, rz(i), n(i)*rz(i+1));
%         ind = maxvol2(crz.');
%         phizz{i} = crz(:,ind);
%         nrmz = 1./min(svd(phizz{i}));
%         phizz{i} = phizz{i}.*nrmz;
        if (~isempty(ifun))
            Jz{i} = indexmerge((1:n(i))', Jz{i+1});
            Jz{i} = Jz{i}(ind,:);
        end;
        if (~isempty(ffun))
            % Elem function is present
            % Interface matrices for TT-fun
            for j=1:nx
                phizx{i,j} = reshape(crX{i,j}, rx(i,j)*n(i), rx(i+1,j));
                phizx{i,j} = phizx{i,j}*phizx{i+1,j};
                phizx{i,j} = reshape(phizx{i,j}, rx(i,j), n(i)*rz(i+1));
                phizx{i,j} = phizx{i,j}(:, ind);
            end;
        end;
        % Extract z indices from y.
        cry = reshape(y{i}, ry(i)*n(i), ry(i+1));
        cry = cry*phizy{i+1};
        cry = reshape(cry, ry(i), n(i)*rz(i+1));
        phizy{i} = cry(:, ind);
        phizy{i} = phizy{i}.*nrms(i);
    end;
end;


if (verb>2)
    Y = cell(d,nswp);
end;

last_sweep = false;
ievalcnt = 0;
fevalcnt = 0;
b = []; % A block size will be stored here
max_dx = 0;
max_dx_prev = inf;
err_raise_cnt = 0;
swp = 1;
dir = 1;
i = 1;
while (swp<=nswp)    
    if (~isempty(ifun))
        % Index function is specified
        % Generate a full index = [Jl, (1:n)', Jr]
        % Jyl to lowest index           n(i) to the middle one,                      Jyr to the senior index
        J = indexmerge(Jy{i}, (1:n(i))', Jy{i+1});
        ievalcnt = ievalcnt + size(J,1);
        if (isempty(b))
            % Take one sample to measure its size =)
            cry1 = ifun(J(1,:));
            b = numel(cry1);
            % Init y{1} by zeros
            y{1} = zeros(ry(1), n(1), ry(2), b);
            ievalcnt = ievalcnt+1;
        end;
        % Actually compute the new core
        if (vec)
            cry = ifun(J);
        else
            % User function is not vectorized - run a loop
            cry = zeros(ry(i)*n(i)*ry(i+1), b);
            for j=1:ry(i)*n(i)*ry(i+1)
                cry(j,:) = ifun(J(j,:));
            end;
        end;
    end;
    if (~isempty(ffun))
        % Compute the X at Y indices
        fx = zeros(ry(i)*n(i)*ry(i+1), nx);
        for j=1:nx
            cr = reshape(crX{i,j}, rx(i,j), n(i)*rx(i+1,j));
            cr = phiyx{i,j}*cr;
            cr = reshape(cr, ry(i)*n(i), rx(i+1,j));
            cr = cr*phiyx{i+1,j};
            fx(:,j) = cr(:);
        end;
        % Call the function
        fevalcnt = fevalcnt+ry(i)*n(i)*ry(i+1);
        fy = ffun(fx); % sizes: rnr x nx -> rnr x d2
        if (isempty(ifun))
            cry = fy;
        else
            cry = cry+fy;
        end;
        if (isempty(b))
            b = size(cry, 2);
            y{1} = zeros(ry(1), n(1), ry(2), b);
        end;
    end;
    
    if (norm(cry, 'fro')==0)
        error('The initial solution is exactly zero. Try a better initial guess\n');
    end;
    
    if (strcmp(trunc_method, 'svd'))
        % Apply interpolation matrices
        cry = reshape(cry, ry(i), n(i)*ry(i+1)*b);
        cry = phiyy{i} \ cry;
        cry = reshape(cry, ry(i)*n(i)*ry(i+1), b);
        cry = cry.';
        cry = reshape(cry, b*ry(i)*n(i), ry(i+1));
        cry = cry/phiyy{i+1};
        cry = reshape(cry, b, ry(i)*n(i)*ry(i+1));
        cry = cry.';
        
        nrms(i) = norm(cry, 'fro');
        cry = cry./nrms(i);
    else
        nrms(i) = 1; % we don't need this in index representation
    end;
    
    y_prev = reshape(y{i}, ry(i)*n(i)*ry(i+1), b);
%     if (strcmp(trunc_method, 'svd'))
%         y_prev = reshape(y_prev, ry(i)*n(i)*ry(i+1), b);
%     else
%         y_prev = reshape(y_prev, ry(i), n(i)*ry(i+1)*b);
%         y_prev = phiyy{i}*y_prev; %
%         y_prev = reshape(y_prev, ry(i)*n(i)*ry(i+1), b);
%         y_prev = y_prev.';
%         y_prev = reshape(y_prev, b*ry(i)*n(i), ry(i+1));
%         y_prev = y_prev*phiyy{i+1};
%         y_prev = reshape(y_prev, b, ry(i)*n(i)*ry(i+1));
%         y_prev = y_prev.';
%         y_prev = y_prev./norm(y_prev, 'fro');
%     end;
    
    dx = norm(cry-y_prev, 'fro'); % /norm(cry,'fro');
    max_dx = max(max_dx,dx);
    
    % Truncation, etc
    if (dir>0)&&(i<d)
        cry = reshape(cry, ry(i)*n(i), ry(i+1)*b);
        if (strcmp(trunc_method, 'svd'))
            [u,s,v]=svd(cry, 'econ');
            s = diag(s);
            r = my_chop2(s, norm(s)*tol/sqrt(d));
            u = u(:,1:r);
            v = diag(s(1:r))*v(:,1:r)';
        else
            % Full-pivot cross should be more accurate
            [u,v]=localcross(cry, tol/sqrt(d));
%             minsz = min(ry(i)*n(i), ry(i+1)*b);
%             u = zeros(ry(i)*n(i), minsz);
%             v = zeros(minsz, ry(i+1)*b);
%             res = cry;
%             val_max = 0;
%             for r=1:minsz
%                 res = reshape(res, [], 1);
%                 [val,piv]=max(abs(res));
%                 val_max = max(val_max, abs(cry(piv)));
%                 piv = tt_ind2sub([ry(i)*n(i), ry(i+1)*b], piv);
%                 if (val<tol*val_max/sqrt(d))
%                     break;
%                 end;
%                 res = reshape(res, ry(i)*n(i), ry(i+1)*b);
%                 u(:,r) = res(:, piv(2));
%                 v(r,:) = res(piv(1), :)/res(piv(1),piv(2));
%                 res = res - u(:,r)*v(r,:);
%             end;
%             u = u(:, 1:r);
%             u = reshape(u, ry(i), n(i)*r); % u
% %             u = phiyy{i}\u;
%             u = reshape(u, ry(i)*n(i), r);
%             v = v(1:r, :);
%             v = reshape(v, r*ry(i+1), b);
%             v = v.';
%             v = reshape(v, b*r, ry(i+1));
% %             v = v/phiyy{i+1};
%             v = reshape(v, b, r*ry(i+1));
%             v = v.';
%             v = reshape(v, r, ry(i+1)*b);
%             % qr u, in case we don't have enrichment
%             [u,rv]=qr(u,0);
%             v = rv*v;
        end;

        if (zrank>0)&&(~last_sweep)
            % Project onto Ez
            crys = u*v;
            crys = reshape(crys, ry(i)*n(i)*ry(i+1), b);
            crys = crys.';
            crys = reshape(crys, b*ry(i)*n(i), ry(i+1));
            crys = crys*phizy{i+1};
            crys = reshape(crys, b, ry(i)*n(i)*rz(i+1));
            crys = crys.';
            cryz = reshape(crys, ry(i), n(i)*rz(i+1)*b);
            cryz = phizy{i}*cryz;
            cryz = reshape(cryz, rz(i)*n(i)*rz(i+1), b);
            
            % Compute Z
            if (~isempty(ifun))
                J = indexmerge(Jz{i}, (1:n(i))', Jz{i+1});
                ievalcnt = ievalcnt + size(J,1);
                if (vec)
                    crz = ifun(J);
                else
                    % User function is not vectorized - run a loop
                    crz = zeros(rz(i)*n(i)*rz(i+1), b);
                    for j=1:rz(i)*n(i)*rz(i+1)
                        crz(j,:) = ifun(J(j,:));
                    end;
                end;
            end;
            if (~isempty(ffun))
                % Compute the X at Y indices
                fx = zeros(rz(i)*n(i)*rz(i+1), nx);
                for j=1:nx
                    cr = reshape(crX{i,j}, rx(i,j), n(i)*rx(i+1,j));
                    cr = phizx{i,j}*cr;
                    cr = reshape(cr, rz(i)*n(i), rx(i+1,j));
                    cr = cr*phizx{i+1,j};
                    fx(:,j) = cr(:);
                end;
                % Call the function
                fevalcnt = fevalcnt+rz(i)*n(i)*rz(i+1);
                fz = ffun(fx); % sizes: rnr x nx -> rnr x d2
                if (isempty(ifun))
                    crz = fz;
                else
                    crz = crz+fz;
                end;
            end;
            crz = crz/nrms(i) - cryz;
%             fprintf('\tswp=%d,i=%d,|y|=%3.3e,|z|=%3.3e,',swp,i,norm(cryz,'fro'),norm(crz,'fro'));
            % Apply interpolation matrices
            crz = reshape(crz, rz(i), n(i)*rz(i+1)*b);
%             crz = phizz{i} \ crz;
            crz = reshape(crz, rz(i)*n(i)*rz(i+1), b);
            crz = crz.';
            crz = reshape(crz, b*rz(i)*n(i), rz(i+1));
%             crz = crz/phizz{i+1};
            crz = reshape(crz, b, rz(i)*n(i)*rz(i+1));
            crz = crz.';
            crz = reshape(crz, rz(i)*n(i), rz(i+1)*b);
            [crz,~]=localcross(crz,0);
%             [crz,~,~]=svd(crz, 'econ');
            crz = crz(:,1:min(size(crz,2),zrank));
            if (zrank2>0)
                crs = [crz, randn(rz(i)*n(i), zrank2)];
                [crz,~]=qr(crs,0);
            end;
            
            % Compute S
            if (~isempty(ifun))
                J = indexmerge(Jy{i}, (1:n(i))', Jz{i+1});
                ievalcnt = ievalcnt + size(J,1);
                if (vec)
                    crs = ifun(J);
                else
                    % User function is not vectorized - run a loop
                    crs = zeros(ry(i)*n(i)*rz(i+1), b);
                    for j=1:ry(i)*n(i)*rz(i+1)
                        crs(j,:) = ifun(J(j,:));
                    end;
                end;
            end;
            if (~isempty(ffun))
                % Compute the X at Y indices
                fx = zeros(ry(i)*n(i)*rz(i+1), nx);
                for j=1:nx
                    cr = reshape(crX{i,j}, rx(i,j), n(i)*rx(i+1,j));
                    cr = phiyx{i,j}*cr;
                    cr = reshape(cr, ry(i)*n(i), rx(i+1,j));
                    cr = cr*phizx{i+1,j};
                    fx(:,j) = cr(:);
                end;
                % Call the function
                fevalcnt = fevalcnt+ry(i)*n(i)*rz(i+1);
                fz = ffun(fx); % sizes: rnr x nx -> rnr x d2
                if (isempty(ifun))
                    crs = fz;
                else
                    crs = crs+fz;
                end;
            end;
            % Apply interpolation matrices - here crys is already in the
            % core variables from the left
            crs = reshape(crs, ry(i), n(i)*rz(i+1)*b);
%             crs = phiyy{i}\crs;
            crs = reshape(crs, ry(i)*n(i)*rz(i+1), b);
            
            crs = crs/nrms(i) - crys;
%             fprintf('|s|=%3.3e\n',norm(crs,'fro'));
            
            crs = crs.';
            crs = reshape(crs, b*ry(i)*n(i), rz(i+1));
%             crs = crs/phizz{i+1};
            crs = reshape(crs, b, ry(i)*n(i)*rz(i+1));
            crs = crs.';
            crs = reshape(crs, ry(i)*n(i), rz(i+1)*b);
            [crs,~]=localcross(crs,0);
%             [crs,~,~]=svd(crs, 'econ');
            crs = crs(:,1:min(size(crs,2),kickrank));
            % Enrichment itself
            [u,rv]=qr([u,crs], 0);
            crs = [v; zeros(size(crs,2), ry(i+1)*b)];
            v = rv*crs;
        end;
        
        if (~strcmp(trunc_method, 'svd'))
            % In the index representation, compute the maxvol immediately
            % and derive u of the form [I ; 0]
            ind = maxvol2(u);
            yy = u(ind,:);
            u = u/yy;
            v = yy*v;
            nrms(i) = 1;
        end;
        
        r = size(u,2);
        v = reshape(v, r*ry(i+1), b);
        v = v.';
        v = reshape(v, b*r, ry(i+1));
        cr2 = reshape(y{i+1}, ry(i+1), n(i+1)*ry(i+2));
        cr2 = v*cr2;
        cr2 = reshape(cr2, b, r*n(i+1)*ry(i+2));
        cr2 = cr2.';
        ry(i+1) = r;        
        y{i} = reshape(u, ry(i), n(i), ry(i+1));
        y{i+1} = reshape(cr2, ry(i+1), n(i+1), ry(i+2), b);
        
        % Maxvols and phis
        if (strcmp(trunc_method, 'svd'))
            u = reshape(u, ry(i), n(i)*ry(i+1));
            u = phiyy{i}*u;
            u = reshape(u, ry(i)*n(i), ry(i+1));
            ind = maxvol2(u);
            phiyy{i+1} = u(ind,:);
            % Extract scales
            nrms(i) = 1./min(svd(phiyy{i+1}));
            phiyy{i+1} = phiyy{i+1}.*nrms(i);
        end;
        if (~isempty(ifun))
            Jy{i+1} = indexmerge(Jy{i}, (1:n(i))');
            Jy{i+1} = Jy{i+1}(ind,:);
        end;
        if (~isempty(ffun))
            for j=1:nx
                phiyx{i+1,j} = reshape(crX{i,j}, rx(i,j), n(i)*rx(i+1,j));
                phiyx{i+1,j} = phiyx{i,j}*phiyx{i+1,j};
                phiyx{i+1,j} = reshape(phiyx{i+1,j}, ry(i)*n(i), rx(i+1,j));
                phiyx{i+1,j} = phiyx{i+1,j}(ind,:);
            end;
        end;
        if (zrank>0)&&(~last_sweep)
            rz(i+1) = size(crz,2);
%             z{i} = crz;
%             crz = reshape(crz, rz(i), n(i)*rz(i+1));
%             crz = phizz{i}*crz;
%             crz = reshape(crz, rz(i)*n(i), rz(i+1));
            ind = maxvol2(crz);
%             phizz{i+1} = crz(ind,:);
%             nrmz = 1./min(svd(phizz{i+1}));
%             phizz{i+1} = phizz{i+1}.*nrmz;
            if (~isempty(ifun))
                Jz{i+1} = indexmerge(Jz{i}, (1:n(i))');
                Jz{i+1} = Jz{i+1}(ind,:);
            end;
            if (~isempty(ffun))
                for j=1:nx
                    phizx{i+1,j} = reshape(crX{i,j}, rx(i,j), n(i)*rx(i+1,j));
                    phizx{i+1,j} = phizx{i,j}*phizx{i+1,j};
                    phizx{i+1,j} = reshape(phizx{i+1,j}, rz(i)*n(i), rx(i+1,j));
                    phizx{i+1,j} = phizx{i+1,j}(ind,:);
                end;
            end;
            phizy{i+1} = reshape(y{i}, ry(i), n(i)*ry(i+1));
            phizy{i+1} = phizy{i}*phizy{i+1};
            phizy{i+1} = reshape(phizy{i+1}, rz(i)*n(i), ry(i+1));
            phizy{i+1} = phizy{i+1}(ind,:);
            phizy{i+1} = phizy{i+1}.*nrms(i);
        end;
    elseif (dir<0)&&(i>1)
        cry = cry.';        
        cry = reshape(cry, b*ry(i), n(i)*ry(i+1));
        if (strcmp(trunc_method, 'svd'))
            [u,s,v]=svd(cry, 'econ');
            s = diag(s);
            r = my_chop2(s, norm(s)*tol/sqrt(d));
            u = u(:,1:r)*diag(s(1:r));
            v = v(:,1:r)';
        else
            % Full-pivot cross should be more accurate
            [v,u]=localcross(cry.', tol/sqrt(d)); % we need orthogonal v
            v = v.';
            u = u.';
%             minsz = min(b*ry(i), n(i)*ry(i+1));
%             u = zeros(b*ry(i), minsz);
%             v = zeros(minsz, n(i)*ry(i+1));
%             res = cry;
%             val_max = 0;
%             for r=1:minsz
%                 res = reshape(res, [], 1);
%                 [val,piv]=max(abs(res));
%                 val_max = max(val_max, abs(cry(piv)));
%                 piv = tt_ind2sub([b*ry(i), n(i)*ry(i+1)], piv);
%                 if (val<tol*val_max/sqrt(d))
%                     break;
%                 end;
%                 res = reshape(res, b*ry(i), n(i)*ry(i+1));
%                 u(:,r) = res(:, piv(2));
%                 v(r,:) = res(piv(1), :)/res(piv(1),piv(2));
%                 res = res - u(:,r)*v(r,:);
%             end;
%             u = u(:, 1:r);
%             u = reshape(u, b, ry(i)*r);
%             u = u.';
%             u = reshape(u, ry(i), r*b);
% %             u = phiyy{i}\u;
%             u = reshape(u, ry(i)*r, b);
%             u = u.';
%             u = reshape(u, b*ry(i), r);
%             v = v(1:r, :);
%             v = reshape(v, r*n(i), ry(i+1));
% %             v = v/phiyy{i+1};
%             v = reshape(v, r, n(i)*ry(i+1));
%             % qr v, in case we don't have enrichment
%             [v,rv]=qr(v.',0);
%             v = v.';
%             u = u*rv.';
        end;
        
        if (zrank>0)&&(~last_sweep)
            % Project onto Ez
            crys = u*v;
            crys = reshape(crys, b, ry(i)*n(i)*ry(i+1));
            crys = crys.';
            crys = reshape(crys, ry(i), n(i)*ry(i+1)*b);
            crys = phizy{i}*crys;
            crys = reshape(crys, rz(i)*n(i)*ry(i+1), b);
            crys = crys.';
            cryz = reshape(crys, b*rz(i)*n(i), ry(i+1));
            cryz = cryz*phizy{i+1};
            cryz = reshape(cryz, b, rz(i)*n(i)*rz(i+1));
            
            % Compute Z
            if (~isempty(ifun))
                J = indexmerge(Jz{i}, (1:n(i))', Jz{i+1});
                ievalcnt = ievalcnt + size(J,1);
                if (vec)
                    crz = fun(J);
                    crz = crz.';
                else
                    % User function is not vectorized - run a loop
                    crz = zeros(b, rz(i)*n(i)*rz(i+1));
                    for j=1:rz(i)*n(i)*rz(i+1)
                        crz(:,j) = fun(J(j,:));
                    end;
                end;
            end;
            if (~isempty(ffun))
                % Compute the X at Y indices
                fx = zeros(rz(i)*n(i)*rz(i+1), nx);
                for j=1:nx
                    cr = reshape(crX{i,j}, rx(i,j), n(i)*rx(i+1,j));
                    cr = phizx{i,j}*cr;
                    cr = reshape(cr, rz(i)*n(i), rx(i+1,j));
                    cr = cr*phizx{i+1,j};
                    fx(:,j) = cr(:);
                end;
                % Call the function
                fevalcnt = fevalcnt+rz(i)*n(i)*rz(i+1);
                fz = ffun(fx); % sizes: rnr x nx -> rnr x d2
                if (isempty(ifun))
                    crz = fz.';
                else
                    crz = crz+fz.';
                end;
            end;
            crz = crz/nrms(i) - cryz;
%             fprintf('\tswp=%d,i=%d,|y|=%3.3e,|z|=%3.3e,',swp,i,norm(cryz,'fro'),norm(crz,'fro'));
            % Apply interpolation matrices
            crz = reshape(crz, b*rz(i)*n(i), rz(i+1));
%             crz = crz/phizz{i+1};            
            crz = reshape(crz, b, rz(i)*n(i)*rz(i+1));
            crz = crz.';
            crz = reshape(crz, rz(i), n(i)*rz(i+1)*b);
%             crz = phizz{i} \ crz;            
            crz = reshape(crz, rz(i)*n(i)*rz(i+1), b);
            crz = crz.';
            crz = reshape(crz, b*rz(i), n(i)*rz(i+1));
            [crz,~]=localcross(crz.',0);
            crz = crz(:,1:min(size(crz,2),zrank));
%             [~,~,crz]=svd(crz, 'econ');
%             crz = conj(crz(:,1:min(size(crz,2),zrank)));
            if (zrank2>0)
                crs = [crz, randn(n(i)*rz(i+1), zrank2)];
                [crz,~]=qr(crs,0);
            end;            
            crz = crz.';
            
            % Compute S
            if (~isempty(ifun))
                J = indexmerge(Jz{i}, (1:n(i))', Jy{i+1});
                ievalcnt = ievalcnt + size(J,1);
                if (vec)
                    crs = fun(J);
                    crs = crs.';
                else
                    % User function is not vectorized - run a loop
                    crs = zeros(b, rz(i)*n(i)*ry(i+1));
                    for j=1:rz(i)*n(i)*ry(i+1)
                        crs(:,j) = fun(J(j,:));
                    end;
                end;
            end;
            if (~isempty(ffun))
                % Compute the X at Y indices
                fx = zeros(rz(i)*n(i)*ry(i+1), nx);
                for j=1:nx
                    cr = reshape(crX{i,j}, rx(i,j), n(i)*rx(i+1,j));
                    cr = phizx{i,j}*cr;
                    cr = reshape(cr, rz(i)*n(i), rx(i+1,j));
                    cr = cr*phiyx{i+1,j};
                    fx(:,j) = cr(:);
                end;
                % Call the function
                fevalcnt = fevalcnt+rz(i)*n(i)*ry(i+1);
                fz = ffun(fx); % sizes: rnr x nx -> rnr x d2
                if (isempty(ifun))
                    crs = fz.';
                else
                    crs = crs+fz.';
                end;
            end;            
            % Apply interpolation matrices
            crs = reshape(crs, b*rz(i)*n(i), ry(i+1));
%             crs = crs/phiyy{i+1};
            crs = reshape(crs, b, rz(i)*n(i)*ry(i+1));
            
            crs = crs/nrms(i) - crys;
%             fprintf('|s|=%3.3e\n', norm(crs,'fro'));

            crs = crs.';
            crs = reshape(crs, rz(i), n(i)*ry(i+1)*b);
%             crs = phizz{i} \ crs;            
            crs = reshape(crs, rz(i)*n(i)*ry(i+1), b);
            crs = crs.';
            crs = reshape(crs, b*rz(i), n(i)*ry(i+1));
            [crs,~]=localcross(crs.', 0);
            crs = crs(:,1:min(size(crs,2),kickrank)).';
%             [~,~,crs]=svd(crs, 'econ');
%             crs = crs(:,1:min(size(crs,2),kickrank))';
            % Enrichment itself
            [v,rv]=qr([v;crs].', 0);
            crs = [u, zeros(b*ry(i), size(crs,1))];
            u = crs*rv.';
            v = v.';
        end;
        
        if (~strcmp(trunc_method, 'svd'))
            ind = maxvol2(v.');
            yy = v(:,ind);
            v = yy\v;
            u = u*yy;
            nrms(i) = 1;
        end;
        
        r = size(v,1);
        u = reshape(u, b, ry(i)*r);
        u = u.';
        u = reshape(u, ry(i), r*b);
        cr2 = reshape(y{i-1}, ry(i-1)*n(i-1), ry(i));
        cr2 = cr2*u;
        ry(i) = r;        
        y{i} = reshape(v, ry(i), n(i), ry(i+1));
        y{i-1} = reshape(cr2, ry(i-1), n(i-1), ry(i), b);
        
        % Maxvols and phis
        if (strcmp(trunc_method, 'svd'))
            v = reshape(v, ry(i)*n(i), ry(i+1));
            v = v*phiyy{i+1};
            v = reshape(v, ry(i), n(i)*ry(i+1));
            ind = maxvol2(v.');
            phiyy{i} = v(:,ind);
            % Extract scales
            nrms(i) = 1./min(svd(phiyy{i}));
            phiyy{i} = phiyy{i}.*nrms(i);
        end;
        if (~isempty(ifun))
            Jy{i} = indexmerge((1:n(i))', Jy{i+1});
            Jy{i} = Jy{i}(ind,:);
        end;
        if (~isempty(ffun))
            % Elem function is present
            % Interface matrices for TT-fun
            for j=1:nx
                phiyx{i,j} = reshape(crX{i,j}, rx(i,j)*n(i), rx(i+1,j));
                phiyx{i,j} = phiyx{i,j}*phiyx{i+1,j};
                phiyx{i,j} = reshape(phiyx{i,j}, rx(i,j), n(i)*ry(i+1));
                phiyx{i,j} = phiyx{i,j}(:, ind);
            end;
        end;
        if (zrank>0)&&(~last_sweep)
            rz(i) = size(crz,1);
%             z{i} = crz;
%             crz = reshape(crz, rz(i)*n(i), rz(i+1));
%             crz = crz*phizz{i+1};
%             crz = reshape(crz, rz(i), n(i)*rz(i+1));
            ind = maxvol2(crz.');
%             phizz{i} = crz(:,ind);
%             nrmz = 1./min(svd(phizz{i}));
%             phizz{i} = phizz{i}.*nrmz;
            if (~isempty(ifun))
                Jz{i} = indexmerge((1:n(i))', Jz{i+1});
                Jz{i} = Jz{i}(ind,:);
            end;
            if (~isempty(ffun))
                % Elem function is present
                % Interface matrices for TT-fun
                for j=1:nx
                    phizx{i,j} = reshape(crX{i,j}, rx(i,j)*n(i), rx(i+1,j));
                    phizx{i,j} = phizx{i,j}*phizx{i+1,j};
                    phizx{i,j} = reshape(phizx{i,j}, rx(i,j), n(i)*rz(i+1));
                    phizx{i,j} = phizx{i,j}(:, ind);
                end;
            end;
            % Extract indices from y.
            cry = reshape(y{i}, ry(i)*n(i), ry(i+1));
            cry = cry*phizy{i+1};
            cry = reshape(cry, ry(i), n(i)*rz(i+1));
            phizy{i} = cry(:, ind);
            phizy{i} = phizy{i}.*nrms(i);   
        end;
    else
%         if (~strcmp(trunc_method, 'svd'))
%             cry = reshape(cry, ry(i), n(i)*ry(i+1)*b);
%             cry = phiyy{i}\cry;
%             cry = reshape(cry, ry(i)*n(i)*ry(i+1), b);
%             cry = cry.';
%             cry = reshape(cry, b*ry(i)*n(i), ry(i+1));
%             cry = cry/phiyy{i+1};
%             cry = reshape(cry, b, ry(i)*n(i)*ry(i+1));
%             cry = cry.';
%         end;
        y{i} = reshape(cry, ry(i), n(i), ry(i+1), b);
    end;
    
    i = i+dir;
    
    % Check the convergence, restart
    if ((i==d+1)||(i==0))
        if (verb>0)
            fprintf('=amen_cross= swp=%d, max_dx=%3.3e, max_rank=%d, #ifun_evals=%d, #ffun_evals=%d\n', swp, max_dx, max(ry), ievalcnt, fevalcnt);
        end;
        
        if (verb>2)
            % Log intermediate solutions
            % Distribute norms equally...
            Y(:,swp) = y;
            nrm1 = exp(sum(log(nrms))/d);
            % ... and plug them into x
            for j=1:d
                Y{j,swp} = Y{j,swp}*nrm1;
            end;
            if (dir>0)
                Y{d,swp} = reshape(Y{d,swp}, ry(d), n(d), b);
            else
                Y{1,swp} = reshape(Y{1,swp}, n(2)*ry(2), b);
                Y{1,swp} = Y{1,swp}.';
                Y{1,swp} = reshape(Y{1,swp}, b, n(2), ry(2));
            end;
            Y{1,swp} = cell2core(tt_tensor,Y(:,swp));
            Y{2,swp} = ievalcnt;
        end;
                      
        if (last_sweep)&&(((exitdir~=0)&&(dir==exitdir))||(exitdir==0))
            break;
        end;
        
        if (max_dx>max_dx_prev)
            err_raise_cnt = err_raise_cnt+1;
        end;

%         if ((max_dx<tol_exit)||(swp==nswp)||(err_raise_cnt>=max_err_jumps))        
        if (max_dx<tol_exit)||(swp==nswp-1)||(err_raise_cnt>=max_err_jumps)
            last_sweep = true;
            % Check if we are going to exit after the wrong sweep.
            % Increase nswp if necessary
            if (exitdir~=0)&&(dir==exitdir)
                nswp = nswp+1;
            end;
            if (err_raise_cnt>=max_err_jumps)
                fprintf('=amen_cross= Will exit since the error stopped decreasing\n');
            end;
%             if (dir==exitdir)
%                 break;
%             else
%                 nswp = nswp+1;
%             end;
        end;
        
        max_dx_prev = max_dx;
        max_dx = 0;
        dir = -dir;
        i = i+dir;
        swp = swp+1;
    end;
end;

% Recover the scales
% Distribute norms equally...
nrms = exp(sum(log(nrms))/d);
% ... and plug them into x
for i=1:d
    y{i} = y{i}*nrms;
end;

if (dir>0)
    y{d} = reshape(y{d}, ry(d), n(d), b);
else
    y{1} = reshape(y{1}, n(1), ry(2), b);
    y{1} = permute(y{1}, [3,1,2]);
end;
y = cell2core(tt_tensor,y);

if (verb>2)
    Y = Y(1:2,1:min(swp,nswp));
end;
end


function [J]=indexmerge(varargin)
% Merges two or three indices in the little-endian manner
sz1 = max(size(varargin{1},1),1);
sz2 = max(size(varargin{2},1),1);
sz3 = 1;
if (nargin>2) % Currently allows only 3
    sz3 = max(size(varargin{3}, 1), 1);
end;
% J1 goes to the fastest index, just copy it
J1 = repmat(varargin{1}, sz2*sz3, 1);
% J2 goes to the middle
J2 = reshape(varargin{2}, 1, []);
J2 = repmat(J2, sz1, 1); % now sz1 ones will be the fastest
J2 = reshape(J2, sz1*sz2, []);
J2 = repmat(J2, sz3, 1);
J = [J1,J2];
if (nargin>2)
    % J3 goes to the slowest
    J3 = reshape(varargin{3}, 1, []);
    J3 = repmat(J3, sz1*sz2, 1); % now sz1 ones will be the fastest
    J3 = reshape(J3, sz1*sz2*sz3, []);
    J = [J,J3];
end;
end

