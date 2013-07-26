function [y]=multifuncrs2(X, funs, eps, varargin)
% Newer Cross approximation of a (vector-)function of several TT-tensors.
%   [Y]=MULTIFUNCRS2(X,FUNS,EPS, VARARGIN)
%   Computes approximation to the functions FUNS(X{1},...,X{N}) with accuracy EPS
%   This version uses a randomized double AMEN enrichment strategy.
%   X should be a cell array of N TT-tensors of equal sizes.
%   The function FUNS should receive a 2d array V of sizes I x N, where the
%   first dimension stays for the reduced set of spatial indices, and  the
%   second is the enumerator of X.
%   The returned sizes should be I x D2, where D2 is the number of
%   components in FUNS. D2 should be either provided as the last (d+1)-th
%   TT-rank of the initial guess, or given explicitly as an option (see
%   below).
%   For example, a linear combination reads FUNS=@(x)(x*W), W is a N x D2
%   matrix (however, you'd better use amen_sum for that...).
%
%   Options are provided in form
%   'PropertyName1',PropertyValue1,'PropertyName2',PropertyValue2 and so
%   on. The parameters are set to default (in brackets in the following)
%   The list of option names and default values are:
%       o y0 - initial approximation [random rank-2 tensor]
%       o nswp - maximal number of DMRG sweeps [10]
%       o rmax - maximal TT rank [Inf]
%       o verb - verbosity level, 0-silent, 1-sweep info, 2-block info [1]
%       o kickrank - the enrichment (AMEn-type) rank for the solution [5]
%       o kickrank2 - the enrichment (rand) rank for the error used above ^^^ [0]
%       o d2 - the last rank of y, that is dim(FUNS) [1]
%       o qr - do (or not) qr before maxvol [false]
%       o eps_exit - stopping accuracy. Try to set it >eps if we stagnate [eps]
%
%********
%   For description of (AMen) adaptive ALS please see
%   Sergey V. Dolgov, Dmitry V. Savostyanov,
%   Alternating minimal energy methods for linear systems in higher dimensions. 
%   Part I: SPD systems, http://arxiv.org/abs/1301.6068,
%   Part II: Faster algorithm and application to nonsymmetric systems, http://arxiv.org/abs/1304.1222
%
%   Use {sergey.v.dolgov, dmitry.savostyanov}@gmail.com for feedback
%*******
%   
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

% Check whether varargin is in fact {varargin}. This happes in restarts
if (numel(varargin)>0)
    testargin = varargin{1};
    if (isa(testargin, 'cell'))
        varargin = testargin;
    end;
end;

nswp = 10;
kickrank = 5;
kickrank2 = 0;
y = [];
verb = 1;
rmax = Inf;
d2 = 1;
do_qr = false;
eps_exit = eps;
restart_it = 0;

for i=1:2:length(varargin)-1
    switch lower(varargin{i})
        case 'nswp'
            nswp=varargin{i+1};
        case 'y0'
            y=varargin{i+1};
        case 'kickrank'
            kickrank=varargin{i+1};
        case 'kickrank2'
            kickrank2=varargin{i+1};            
        case 'rmax'
            rmax=varargin{i+1};            
        case 'verb'
            verb=varargin{i+1};
%         case 'kicktype'
%             kicktype=varargin{i+1};            
%         case 'pcatype'
%             pcatype=varargin{i+1};
%         case 'trunctype'
%             trunctype=varargin{i+1};            
        case 'd2'
            d2=varargin{i+1};            
        case 'qr'
            do_qr = varargin{i+1};
        case 'eps_exit'
            eps_exit = varargin{i+1};            
        case 'restart_it'
            restart_it = varargin{i+1};
        otherwise
            error('Unrecognized option: %s\n',varargin{i});
    end
end

nx = numel(X);
d = X{1}.d;
n = X{1}.n;
rx = zeros(d+1,nx);
crX = cell(d,nx);
for i=1:nx
    rx(:,i) = X{i}.r;
    crX(:,i) = core2cell(X{i});
end;

wasrand = false;
if (isempty(y))
    ry = d2*ones(d+1,1); ry(1)=1;
    y = tt_rand(n, d, ry);
    wasrand = true;
end;
ry = y.r;
cry = core2cell(y);

% error vector
z = tt_rand(n,d,kickrank);
rz = z.r;
z = core2cell(z);

% Interface matrices - for solution
Ry = cell(d+1,1);
Ry{1} = 1; Ry{d+1}=1;
% For error:
% these will be for stuffing z indices into y
Ryz = cell(d+1,1);
Ryz{1} = 1; Ryz{d+1}=1;
% these are for z itself
Rz = cell(d+1,1);
Rz{1} = 1; Rz{d+1}=1;
Rx = cell(d+1,nx);
% these are for stuffing z indices into x
Rxz = cell(d+1,nx);
for i=1:nx
    Rx{1,i}=1; Rx{d+1,i}=1;
    Rxz{1,i}=1; Rxz{d+1,i}=1;
end;

block_order = [+(d), -(d)];

% Orth
for i=1:d-1
    cr = cry{i}; % r1,n,d2,r2
    cr = reshape(cr, ry(i)*n(i), ry(i+1));
    [cr, rv]=qr(cr, 0);    
    cr2 = cry{i+1};
    cr2 = reshape(cr2, ry(i+1), n(i+1)*ry(i+2));
    cr2 = rv*cr2;
    ry(i+1) = size(cr, 2);
    cr = reshape(cr, ry(i), n(i), ry(i+1));
    cry{i+1} = reshape(cr2, ry(i+1), n(i+1), ry(i+2));
    cry{i} = cr;

    % Interfaces for solution
    % Interface matrix for Y        
    Ry{i+1} = Ry{i}*reshape(cr, ry(i), n(i)*ry(i+1));
    Ry{i+1} = reshape(Ry{i+1}, ry(i)*n(i), ry(i+1));
    if (wasrand)
        curind = [];
        while numel(curind)<ry(i+1)
            curind = [curind; ceil(rand(ry(i+1), 1)*(n(i)*ry(i)))];
            curind = unique(curind);
        end;
        curind = curind(1:ry(i+1));
    else
        curind = maxvol2(Ry{i+1},'qr',do_qr);
    end;    
    Ry{i+1} = Ry{i+1}(curind, :);
    % Interface matrices for X
    for j=1:nx
        Rx{i+1,j} = reshape(crX{i,j}, rx(i,j), n(i)*rx(i+1,j));
        Rx{i+1,j} = Rx{i,j}*Rx{i+1,j};
        Rx{i+1,j} = reshape(Rx{i+1,j}, ry(i)*n(i), rx(i+1,j));
        Rx{i+1,j} = Rx{i+1,j}(curind, :);
    end;   
    
    % error for kick
    crz = z{i}; % r1,n,d2,r2
    crz = reshape(crz, rz(i)*n(i), rz(i+1));
    [crz, rv]=qr(crz, 0);    
    cr2 = z{i+1};
    cr2 = reshape(cr2, rz(i+1), n(i+1)*rz(i+2));
    cr2 = rv*cr2;
    rz(i+1) = size(crz, 2);
    crz = reshape(crz, rz(i), n(i), rz(i+1));
    z{i+1} = reshape(cr2, rz(i+1), n(i+1), rz(i+2));
    z{i} = crz;
    
    % Interfaces for error
    Rz{i+1} = Rz{i}*reshape(crz, rz(i), n(i)*rz(i+1));
    Rz{i+1} = reshape(Rz{i+1}, rz(i)*n(i), rz(i+1));
    Ryz{i+1} = Ryz{i}*reshape(cr, ry(i), n(i)*ry(i+1));
    Ryz{i+1} = reshape(Ryz{i+1}, rz(i)*n(i), ry(i+1));    
    % Pick random initial indices
    curind = [];
    while numel(curind)<rz(i+1)
        curind = [curind; ceil(rand(rz(i+1), 1)*(n(i)*rz(i)))];
        curind = unique(curind);
    end;
    curind = curind(1:rz(i+1));
    
    Ryz{i+1} = Ryz{i+1}(curind, :);
    Rz{i+1} = Rz{i+1}(curind, :);
    % Interface matrices for X
    for j=1:nx
        Rxz{i+1,j} = reshape(crX{i,j}, rx(i,j), n(i)*rx(i+1,j));
        Rxz{i+1,j} = Rxz{i,j}*Rxz{i+1,j};
        Rxz{i+1,j} = reshape(Rxz{i+1,j}, rz(i)*n(i), rx(i+1,j));
        Rxz{i+1,j} = Rxz{i+1,j}(curind, :);
    end;    
end;


d2 = ry(d+1);
ry(d+1) = 1;
cry{d} = permute(cry{d}, [3,1,2]); % d2, rd, nd

swp = 1;

max_dy = 0;

cur_order = block_order;
order_index = 2;
i = d;
dir = sign(cur_order(order_index));

% DMRG sweeps
while (swp<=nswp)||(dir>0)
    
    oldy = reshape(cry{i}, d2*ry(i)*n(i)*ry(i+1), 1);
    
    % Compute the X superblocks
    curbl = zeros(ry(i)*n(i)*ry(i+1), nx);
    for j=1:nx
        cr = reshape(crX{i,j}, rx(i,j), n(i)*rx(i+1,j));
        cr = Rx{i,j}*cr;
        cr = reshape(cr, ry(i)*n(i), rx(i+1,j));
        cr = cr*Rx{i+1,j};
        curbl(:,j) = cr(:);
    end;
    % Call the function
    newy = funs(curbl); % sizes: rnr x nx -> rnr x d2
    % Multiply with inverted Ry
    newy = reshape(newy, ry(i), n(i)*ry(i+1)*d2);
    newy = (Ry{i}) \ newy;
    newy = reshape(newy, ry(i)*n(i)*ry(i+1), d2);
    newy = reshape(newy.', d2*ry(i)*n(i), ry(i+1));
    newy = newy / (Ry{i+1});
    newy = reshape(newy, d2*ry(i)*n(i)*ry(i+1), 1);
    
    dy = norm(newy(:)-oldy)/norm(newy(:));
    if (isinf(dy))
        fprintf('Bad initial indices, the solution is exactly zero. Restarting...\n');
        if (restart_it>8)
            error('Restart counts exceeded. Very bad initial guess or function, giving up...\n');
        else
            y = multifuncrs2(X,funs,eps,[varargin, {'restart_it', restart_it+1}]);
            return;
        end;
    end;
    max_dy = max(max_dy, dy);

    % Truncation
    if (dir>0) % left-to-right
        newy = reshape(newy, d2, ry(i)*n(i)*ry(i+1));
        newy = reshape(newy.', ry(i)*n(i), ry(i+1)*d2);
    else
        newy = reshape(newy, d2*ry(i), n(i)*ry(i+1));
    end;

    if (kickrank>=0)
        [u,s,v]=svd(newy, 'econ');
        s = diag(s);
        r = my_chop2(s, eps/sqrt(d)*norm(s));
        r = min(r, rmax);
    else
        if (dir>0)
            [u,v]=qr(newy, 0);
            v=v';
            r = size(u,2);
            s = ones(r,1);
        else
            [v,u]=qr(newy.', 0);
            v=conj(v);
            u=u.';
            r = size(u,2);
            s = ones(r,1);
        end;
    end;

    if (verb>1)
    	fprintf('=multifuncrs=   block %d{%d}, dy: %3.3e, r: %d\n', i, dir, dy, r);
    end;    
    
    % Kicks and interfaces
    if (dir>0)&&(i<d) % left-to-right, kickrank, etc
        u = u(:,1:r);
        v = conj(v(:,1:r))*diag(s(1:r));

        % kick
        radd = 0; rv = 1;
        if (kickrank>0)
            % AMEn kick. See also Cross-3D @Savostyanov, and randomized
            % TT-cross @Zheltkov.
            % However, in this method (AMEn) the residual is stored and
            % approximated separately, which is not exactly (??) Cross-3D.
            
            % Compute the function at residual indices
            curbl_y = zeros(ry(i)*n(i)*rz(i+1), nx);
            curbl_z = zeros(rz(i)*n(i)*rz(i+1), nx);
            for j=1:nx
                % for kick
                cr = reshape(crX{i,j}, rx(i,j), n(i)*rx(i+1,j));
                cr = Rx{i,j}*cr;
                cr = reshape(cr, ry(i)*n(i), rx(i+1,j));
                cr = cr*Rxz{i+1,j};
                curbl_y(:,j) = cr(:);
                % for z update
                cr = reshape(crX{i,j}, rx(i,j), n(i)*rx(i+1,j));
                cr = Rxz{i,j}*cr;
                cr = reshape(cr, rz(i)*n(i), rx(i+1,j));
                cr = cr*Rxz{i+1,j};
                curbl_z(:,j) = cr(:);                
            end;
            % Call the function
            zy = funs(curbl_y); % sizes: rnr x nx -> rnr x d2
            zz = funs(curbl_z); % sizes: rnr x nx -> rnr x d2
            % Assemble y at z indices (sic!) and subtract
            dzy = reshape(u*v.', ry(i)*n(i)*ry(i+1), d2);
            dzy = reshape(dzy.', d2*ry(i)*n(i), ry(i+1));
            % Cast dzy from core items to samples at right indices
            dzy = dzy*Ryz{i+1};
            dzy = reshape(dzy, d2, ry(i)*n(i)*rz(i+1));
            dzy = dzy.';
            % zy still requires casting from samples to core entries
            zy = reshape(zy, ry(i), n(i)*rz(i+1)*d2);
            zy = Ry{i}\zy;
            zy = reshape(zy, ry(i)*n(i)*rz(i+1), d2);
            zy = zy - dzy; % Core elements from the left, indices from right
            % Prepare pure Z:
            % cast dzy from core items to samples at left idx as well
            dzy = reshape(dzy, ry(i), n(i)*rz(i+1)*d2);
            dzy = Ryz{i}*dzy;
            dzy = reshape(dzy, rz(i)*n(i)*rz(i+1), d2);
            zz = zz - dzy; % Samples from both sides

            % Interp all remaining samples into core elements
            % ...for kick
            zy = reshape(zy.', d2*ry(i)*n(i), rz(i+1));
            zy = zy / Rz{i+1};
            zy = reshape(zy, d2, ry(i)*n(i)*rz(i+1));
            zy = reshape(zy.', ry(i)*n(i), rz(i+1)*d2);
            % SVD to eliminate d2 and possibly overestimated rz
            [zy,sz,vz]=svd(zy, 'econ');
            zy = zy(:, 1:min(kickrank,size(zy,2)));
            % ...for z update
            zz = reshape(zz, rz(i), n(i)*rz(i+1)*d2);
            zz = Rz{i} \ zz;
            zz = reshape(zz, rz(i)*n(i)*rz(i+1), d2);
            zz = reshape(zz.', d2*rz(i)*n(i), rz(i+1));
            zz = zz / Rz{i+1};
            zz = reshape(zz, d2,rz(i)*n(i)*rz(i+1));
            zz = reshape(zz.', rz(i)*n(i), rz(i+1)*d2);
            [zz,sz,vz]=svd(zz, 'econ');            
            zz = zz(:, 1:min(kickrank,size(zz,2)));
            % second random kick - for z            
            zz = [zz, randn(rz(i)*n(i), kickrank2)];
            
            [u,rv]=qr([u,zy], 0);
            radd = size(zy,2);
        end;
        v = [v, zeros(ry(i+1)*d2, radd)];
        v = rv*v.';
        r = size(u,2);

        cr2 = cry{i+1};
        cr2 = reshape(cr2, ry(i+1), n(i+1)*ry(i+2));
        v = reshape(v, r*ry(i+1), d2);
        v = reshape(v.', d2*r, ry(i+1));
        v = v*cr2; % size r+radd, n2, r3

        ry(i+1) = r;

        u = reshape(u, ry(i), n(i), r);
        v = reshape(v, d2, r, n(i+1), ry(i+2));

        % Stuff back
        cry{i} = u;
        cry{i+1} = v;
        
        % Update kick
        [zz,rv]=qr(zz,0);        
        rz(i+1) = size(zz,2);
        z{i} = reshape(zz, rz(i), n(i), rz(i+1));
        % z{i+1} is recomputed from scratch, we don't need it now        
        
        % Recompute left interface matrices
        % Interface matrix for Y
        Ry{i+1} = Ry{i}*reshape(u, ry(i), n(i)*ry(i+1));
        Ry{i+1} = reshape(Ry{i+1}, ry(i)*n(i), ry(i+1));
        curind = maxvol2(Ry{i+1},'qr',do_qr);
        Ry{i+1} = Ry{i+1}(curind, :);
        % Interface matrices for X
        for j=1:nx
            Rx{i+1,j} = reshape(crX{i,j}, rx(i,j), n(i)*rx(i+1,j));
            Rx{i+1,j} = Rx{i,j}*Rx{i+1,j};
            Rx{i+1,j} = reshape(Rx{i+1,j}, ry(i)*n(i), rx(i+1,j));
            Rx{i+1,j} = Rx{i+1,j}(curind, :);
        end;
        % for kick
        Ryz{i+1} = Ryz{i}*reshape(u, ry(i), n(i)*ry(i+1));
        Ryz{i+1} = reshape(Ryz{i+1}, rz(i)*n(i), ry(i+1));
        Rz{i+1} = Rz{i}*reshape(zz, rz(i), n(i)*rz(i+1));
        Rz{i+1} = reshape(Rz{i+1}, rz(i)*n(i), rz(i+1));        
        curind = maxvol2(Rz{i+1},'qr',do_qr);
        Ryz{i+1} = Ryz{i+1}(curind, :);
        Rz{i+1} = Rz{i+1}(curind, :);
        % Interface matrices for X
        for j=1:nx
            Rxz{i+1,j} = reshape(crX{i,j}, rx(i,j), n(i)*rx(i+1,j));
            Rxz{i+1,j} = Rxz{i,j}*Rxz{i+1,j};
            Rxz{i+1,j} = reshape(Rxz{i+1,j}, rz(i)*n(i), rx(i+1,j));
            Rxz{i+1,j} = Rxz{i+1,j}(curind, :);
        end;        
    elseif (dir<0)&&(i>1) % right-to-left
        u = u(:,1:r)*diag(s(1:r));
        v = conj(v(:,1:r));
        % kick
        radd = 0; rv = 1;
        if (kickrank>0)
            % AMEn kick
            % Compute the function at residual indices
            curbl_y = zeros(rz(i)*n(i)*ry(i+1), nx);
            curbl_z = zeros(rz(i)*n(i)*rz(i+1), nx);
            for j=1:nx
                % for kick
                cr = reshape(crX{i,j}, rx(i,j), n(i)*rx(i+1,j));
                cr = Rxz{i,j}*cr;
                cr = reshape(cr, rz(i)*n(i), rx(i+1,j));
                cr = cr*Rx{i+1,j};
                curbl_y(:,j) = cr(:);
                % for z update
                cr = reshape(crX{i,j}, rx(i,j), n(i)*rx(i+1,j));
                cr = Rxz{i,j}*cr;
                cr = reshape(cr, rz(i)*n(i), rx(i+1,j));
                cr = cr*Rxz{i+1,j};
                curbl_z(:,j) = cr(:);                
            end;
            % Call the function
            zy = funs(curbl_y); % sizes: rnr x nx -> rnr x d2
            zz = funs(curbl_z); % sizes: rnr x nx -> rnr x d2
            % Assemble y at z indices (sic!) and subtract
            dzy = reshape(u*v.', ry(i), n(i)*ry(i+1)*d2);
            dzy = Ryz{i}*dzy;
            dzy = reshape(dzy, rz(i)*n(i)*ry(i+1), d2);
            % zy still requires casting from samples to core entries
            zy = zy.';
            zy = reshape(zy, d2*rz(i)*n(i), ry(i+1));
            zy = zy/Ry{i+1};
            zy = reshape(zy, d2, rz(i)*n(i)*ry(i+1));
            zy = zy.';
            zy = zy - dzy; % samples at left, core elems at right
            dzy = reshape(dzy.', d2*rz(i)*n(i), ry(i+1));
            dzy = dzy*Ryz{i+1};
            dzy = reshape(dzy, d2, rz(i)*n(i)*rz(i+1));
            zz = zz - (dzy.'); % samles everywhere
            
            % Cast sample indices to core elements
            % ...for kick
            zy = reshape(zy, rz(i), n(i)*ry(i+1)*d2);
            zy = (Rz{i}) \ zy;
            zy = reshape(zy, rz(i)*n(i)*ry(i+1), d2);
            zy = zy.';
            zy = reshape(zy, d2*rz(i), n(i)*ry(i+1));
            [zu,zs,zy]=svd(zy, 'econ');
            zy = conj(zy(:, 1:min(kickrank,size(zy,2))));
            % ...for z update
            zz = reshape(zz, rz(i), n(i)*rz(i+1)*d2);
            zz = (Rz{i}) \ zz;
            zz = reshape(zz, rz(i)*n(i)*rz(i+1), d2);
            zz = reshape(zz.', d2*rz(i)*n(i), rz(i+1));
            zz = zz / (Rz{i+1});
            zz = reshape(zz, d2*rz(i), n(i)*rz(i+1));
            [zu,zs,zz]=svd(zz, 'econ');
            zz = conj(zz(:,1:min(kickrank, size(zz,2))));
            zz = [zz, randn(n(i)*rz(i+1), kickrank2)];
            
            [v,rv]=qr([v,zy], 0);
            radd = size(zy,2);
        end;
        u = [u, zeros(d2*ry(i), radd)];
        u = u*rv.';
        r = size(v,2);
        
        cr2 = cry{i-1};
        cr2 = reshape(cr2, ry(i-1)*n(i-1), ry(i));
        u = reshape(u, d2, ry(i)*r);
        u = reshape(u.', ry(i), r*d2);
        u = cr2*u;
               
        u = reshape(u, ry(i-1)*n(i-1)*r, d2);
        u = reshape(u.', d2, ry(i-1), n(i-1), r);
        v = reshape(v.', r, n(i), ry(i+1));
                
        % Stuff back
        ry(i) = r;
        cry{i-1} = u;
        cry{i} = v;
        
        % kick
        [zz,rv]=qr(zz,0);
        rz(i) = size(zz,2);
        zz = reshape(zz.', rz(i), n(i), rz(i+1));
        z{i} = zz;
        % z{i-1} is recomputed from scratch, we don't need it now
        
        % Recompute left interface matrices
        % Interface matrix for Y
        Ry{i} = reshape(v, ry(i)*n(i), ry(i+1))*Ry{i+1};
        Ry{i} = reshape(Ry{i}, ry(i), n(i)*ry(i+1));
        curind = maxvol2(Ry{i}.','qr',do_qr);
        Ry{i} = Ry{i}(:, curind);
        % Interface matrices for X
        for j=1:nx
            Rx{i,j} = reshape(crX{i,j}, rx(i,j)*n(i), rx(i+1,j));
            Rx{i,j} = Rx{i,j}*Rx{i+1,j};
            Rx{i,j} = reshape(Rx{i,j}, rx(i,j), n(i)*ry(i+1));
            Rx{i,j} = Rx{i,j}(:, curind);
        end;
        % for kick
        Rz{i} = reshape(zz, rz(i)*n(i), rz(i+1))*Rz{i+1};
        Rz{i} = reshape(Rz{i}, rz(i), n(i)*rz(i+1));
        Ryz{i} = reshape(v, ry(i)*n(i), ry(i+1))*Ryz{i+1};
        Ryz{i} = reshape(Ryz{i}, ry(i), n(i)*rz(i+1));        
        curind = maxvol2(Rz{i}.','qr',do_qr);
        Ryz{i} = Ryz{i}(:, curind);
        Rz{i} = Rz{i}(:, curind);
        % Interface matrices for X
        for j=1:nx
            Rxz{i,j} = reshape(crX{i,j}, rx(i,j)*n(i), rx(i+1,j));
            Rxz{i,j} = Rxz{i,j}*Rxz{i+1,j};
            Rxz{i,j} = reshape(Rxz{i,j}, rx(i,j), n(i)*rz(i+1));
            Rxz{i,j} = Rxz{i,j}(:, curind);
        end;        
    elseif ((dir>0)&&(i==d))
        % Just stuff back the last core
        newy = u(:,1:r)*diag(s(1:r))*v(:,1:r)';
        newy = reshape(newy, ry(i)*n(i)*ry(i+1), d2);
        cry{i} = reshape(newy.', d2, ry(i), n(i), ry(i+1));
    elseif ((dir<0)&&(i==1))
        % Just stuff back the last core
        newy = u(:,1:r)*diag(s(1:r))*v(:,1:r)';
        newy = reshape(newy, d2, ry(i), n(i), ry(i+1));
        cry{i} = newy;        
    end;
    
    
    i = i+dir;
    % Reversing, residue check, etc
    cur_order(order_index) = cur_order(order_index) - dir;
    % New direction
    if (cur_order(order_index)==0)
        order_index = order_index+1;

        if (verb>0)
            fprintf('=multifuncrs= sweep %d{%d}, max_dy: %3.3e, erank: %g\n', swp, order_index-1, max_dy, sqrt(ry(1:d)'*(n.*ry(2:d+1))/sum(n)));
        end;

        if (max_dy<eps_exit)&&(dir>0)
            break;
        end;

        if (order_index>numel(cur_order)) % New global sweep
            cur_order = block_order;
            order_index = 1;

            max_dy = 0;
            swp = swp+1;
        end;

        dir = sign(cur_order(order_index));
        i = i+dir;
    end;
end

cry{d} = permute(cry{d}, [2,3,1]); % d2 is r(d+1)
y = cell2core(y, cry);

end
