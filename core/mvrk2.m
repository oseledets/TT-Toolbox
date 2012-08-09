function [y]=mvrk2(a,x,tol,varargin)
%ALS-WDB method for the MatVec approximation in QTT-Tucker format
%   [Y]=MVRK2(A,X,TOL,VARARGIN) Attempts to approximate the product
%   Y = A*X with accuracy EPS using the AMR (Wedderburn) iteration.
%   Matrix A has to be given in the QTT-Tucker, vector X should be
%   given in the QTT-Tucker format also. Options are provided in form
%   'PropertyName1',PropertyValue1,'PropertyName2',PropertyValue2 and so
%   on. The parameters are set to default (in brackets in the following) 
%   The list of option names and default values are:
%       o y0 - initial approximation [tensor of all ones] 
%       o nswp - maximal number of AMR sweeps [10]
%       o verb - verbosity level, 0-silent, 1-sweep info, 2-block info [1]
%       o kickrank - rank-increasing parameter [5]
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

% Inner parameters

nswp=10;
als_tol_low = 2;
als_tol_high = 10;
als_iters = 2;

% rmax=1000;
verb=1;
kickrank = 5;
y=[];

for i=1:2:length(varargin)-1
    switch lower(varargin{i})
        case 'nswp'
            nswp=varargin{i+1};
%         case 'rmax'
%             rmax=varargin{i+1};
        case 'y0'
            y=varargin{i+1};
        case 'verb'
            verb=varargin{i+1};
        case 'als_tol_high'
            als_tol_high=varargin{i+1};                        
        case 'als_tol_low'
            als_tol_low=varargin{i+1};
        case 'als_iters'
            als_iters=varargin{i+1};            
        case 'kickrank'
            kickrank=varargin{i+1};
            
        otherwise
            error('Unrecognized option: %s\n',varargin{i});
    end
end

tol2 = tol;

d = x.core.d;
xc = core2cell(x.core);
rxc = x.core.r;
m = cell(d,1);
L = zeros(d,1);
xf = cell(d,1);
rxf = cell(d,1);
for i=1:d
    m{i} = x.tuck{i}.n;
    L(i) = x.tuck{i}.d;
    rxf{i} = x.tuck{i}.r;
    rxf{i}(L(i)+2) = rxc(i+1);
    xf{i} = cell(L(i)+1,1);
    xf{i}(1:L(i)) = core2cell(x.tuck{i});
    xf{i}{L(i)+1} = xc{i};
end;

ac = core2cell(a.core);
rac = a.core.r;
n = cell(d,1);
af = cell(d,1);
raf = cell(d,1);
for i=1:d
    raf{i} = a.tuck{i}.r;
    n{i} = a.tuck{i}.n;
    raf{i}(L(i)+2) = rac(i+1);    
    af{i} = cell(L(i)+1,1);
    af{i}(1:L(i)) = core2cell(a.tuck{i});
    af{i}{L(i)+1} = ac{i};
end;

yf = cell(d,1);
ryf = cell(d,1);
if (isempty(y))
    y = qtt_tucker;
    yc = tt_ones(1, d);
    y.core = yc;
    ryc = yc.r;
    yc = core2cell(yc);
    y.tuck = cell(d,1);
    for i=1:d
        ryf{i} = ones(L(i)+2,1);
        n{i}(L(i)+1) = 1;
        yf{i} = cell(L(i)+1,1);
        y.tuck{i} = tt_ones(n{i}(1:L(i)));
        yf{i}(1:L(i)) = core2cell(y.tuck{i});
        yf{i}{L(i)+1}=yc{i};
    end;
else
    yc = core2cell(y.core);
    ryc = y.core.r;
    for i=1:d
        ryf{i} = y.tuck{i}.r;
        ryf{i}(L(i)+2)=ryc(i+1);
        n{i}(L(i)+1) = ryc(i);
        yf{i} = cell(L(i)+1,1);
        yf{i}(1:L(i)) = core2cell(y.tuck{i});
        yf{i}{L(i)+1}=yc{i};
    end;    
end;

max_dx = 0;
max_dx_prev=Inf;
regurg_cnt=0;

phiaf = cell(d,1);
for i=1:d
    phiaf{i} = cell(L(i)+2,1);
    phiaf{i}{1}=1;
end;
phiac = cell(d+1,1); phiac{1}=1; phiac{d+1}=1;

acp = cell(d,1);
% xcp = cell(d,1);

last_sweep = false;
max_frank = 0;

for swp=1:nswp
    
    % Orthog + phi, Factors
    for i=1:d
        % Permute core blocks: rc1 is mode index now
        yf{i}{L(i)+1} = permute(yf{i}{L(i)+1}, [2,1,3]); % rxt,rxc1,rxc2
        for j=1:L(i)
            cr = yf{i}{j};
            cr = reshape(cr, ryf{i}(j)*n{i}(j), ryf{i}(j+1));
            [cr,rv]=qr(cr,0);
            cr2 = yf{i}{j+1};
            cr2 = reshape(cr2, ryf{i}(j+1), n{i}(j+1)*ryf{i}(j+2));
            cr2 = rv*cr2;
            ryf{i}(j+1) = size(cr, 2);
            cr = reshape(cr, ryf{i}(j), n{i}(j), ryf{i}(j+1));
            yf{i}{j} = cr;
            yf{i}{j+1} = reshape(cr2, ryf{i}(j+1), n{i}(j+1), ryf{i}(j+2));
            
            phiaf{i}{j+1} = compute_next_Phi(phiaf{i}{j}, cr, af{i}{j}, xf{i}{j}, 'lr');
        end;
        % core-projected matrix block
        acp{i} = permute(phiaf{i}{L(i)+1}, [2, 1, 3]);
        acp{i} = reshape(acp{i}, raf{i}(L(i)+1), ryf{i}(L(i)+1)*rxf{i}(L(i)+1));
        acp{i} = reshape(permute(af{i}{L(i)+1}, [1,3,2]), rac(i)*rac(i+1), raf{i}(L(i)+1)) * acp{i};
        acp{i} = reshape(acp{i}, rac(i), rac(i+1), ryf{i}(L(i)+1), rxf{i}(L(i)+1));
        acp{i} = permute(acp{i}, [1, 3, 4, 2]);
%         % core-projected rhs block
%         ycp{i} = phiyf{i}{L(i)+1}.';
%         ycp{i} = reshape(permute(yf{i}{L(i)+1}, [1,3,2]), ryc(i)*ryc(i+1), ryf{i}(L(i)+1)) * ycp{i};
%         ycp{i} = reshape(ycp{i}, ryc(i), ryc(i+1), rxf{i}(L(i)+1));
%         ycp{i} = permute(ycp{i}, [1, 3, 2]);
        % return xf{i}{L{i}+1} to the core state
        yf{i}{L(i)+1} = permute(yf{i}{L(i)+1}, [2,1,3]); % rxc1,rxt,rxc2
    end;
    
    % ALS over core
    nc = zeros(d,1);
    mc = zeros(d,1);
    for i=1:d
        nc(i)=ryf{i}(L(i)+1);
        mc(i)=rxf{i}(L(i)+1);
    end;
    % Orthog
    for i=1:d-1
        cr = yf{i}{L(i)+1};
        cr = reshape(cr, ryc(i)*nc(i), ryc(i+1));
        [cr,rv]=qr(cr,0);
        cr2 = yf{i+1}{L(i+1)+1};
        cr2 = reshape(cr2, ryc(i+1), nc(i+1)*ryc(i+2));
        cr2 = rv*cr2;
        ryc(i+1) = size(cr,2);
        % Update Factor-sizes too
        n{i+1}(L(i+1)+1)=ryc(i+1);
        ryf{i}(L(i)+2)=ryc(i+1);
        cr = reshape(cr, ryc(i), nc(i), ryc(i+1));
        yf{i}{L(i)+1} = cr;
        yf{i+1}{L(i+1)+1} = reshape(cr2, ryc(i+1), nc(i+1), ryc(i+2));
        
        phiac{i+1} = compute_next_Phi(phiac{i}, cr, acp{i}, xf{i}{L(i)+1}, 'lr');
    end;
    % ALS
    for i=d:-1:1
        Phi1 = phiac{i};
        A1 = acp{i};
        x1 = xf{i}{L(i)+1};
        Phi2 = phiac{i+1};
        
        y_prev = reshape(yf{i}{L(i)+1}, ryc(i)*nc(i)*ryc(i+1), 1);
        
        dir = -1;
        if (i==1); dir=0; end;
        local_kickrank = kickrank;
        if (last_sweep); local_kickrank=0; end;
        [u,v,max_dx]=local_proj(Phi1,A1,Phi2, x1, tol2/sqrt(sum(L+1)), y_prev,  max_dx, dir, local_kickrank, verb);
        
        if (i>1)
            cr2 = yf{i-1}{L(i-1)+1};
            cr2 = reshape(cr2, ryc(i-1)*nc(i-1), ryc(i));
            cr2 = cr2*u;
            ryc(i) = size(v,2);
            % Update Factor-sizes
            n{i}(L(i)+1)=ryc(i);
            ryf{i-1}(L(i-1)+2)=ryc(i);
            v = reshape(v.', ryc(i), nc(i), ryc(i+1));
            yf{i}{L(i)+1} = v;
            yf{i-1}{L(i-1)+1} = reshape(cr2, ryc(i-1), nc(i-1), ryc(i));
            
            phiac{i} = compute_next_Phi(phiac{i+1}, v, acp{i}, x1, 'rl');
        else
            v = u*(v.');
            yf{i}{L(i)+1} = reshape(v, ryc(i), nc(i), ryc(i+1));
        end;
    end;
    
    
    % Optimization over factors
    for i=1:d
        % ALS L->1, the orth. holds
        % Copy the extended factor data to the current array
        cury = yf{i};
        % Reshape the core block to be the factor one
        cury{L(i)+1} = permute(cury{L(i)+1}, [2, 1, 3]);
        
        cura = af{i};
        % last block has to be convolved with phiac-left
        cura{L(i)+1} = permute(phiac{i}, [1,3,2]);
        cura{L(i)+1} = reshape(cura{L(i)+1}, ryc(i)*rxc(i), rac(i));
        cura{L(i)+1} = cura{L(i)+1}*reshape(af{i}{L(i)+1}, rac(i), raf{i}(L(i)+1)*rac(i+1));
        cura{L(i)+1} = reshape(cura{L(i)+1}, ryc(i), rxc(i), raf{i}(L(i)+1), rac(i+1));
        cura{L(i)+1} = permute(cura{L(i)+1}, [3, 1,2, 4]); % tucker rank is now the left
        
        curx = xf{i};
        curx{L(i)+1} = permute(curx{L(i)+1}, [2, 1, 3]); % tucker rank is now the left        
        
        % The L+2-th phi comes from the core
        phiaf{i}{L(i)+2} = phiac{i+1}; % sizes rxc2,ra2,rxc2. rxc2 - our last "rank" index
        
        % First, orthogonality L->1
        for j=(L(i)+1):-1:2
            cr = cury{j};
            cr = reshape(cr, ryf{i}(j), n{i}(j)*ryf{i}(j+1));
            [cr,rv]=qr(cr.',0);
            cr2 = cury{j-1};
            cr2 = reshape(cr2, ryf{i}(j-1)*n{i}(j-1), ryf{i}(j));
            cr2 = cr2*(rv.');
            ryf{i}(j) = size(cr, 2);
            cr = reshape(cr.', ryf{i}(j), n{i}(j), ryf{i}(j+1));
            cury{j} = cr;
            cury{j-1} = reshape(cr2, ryf{i}(j-1), n{i}(j-1), ryf{i}(j));
            
            phiaf{i}{j} = compute_next_Phi(phiaf{i}{j+1}, cr, cura{j}, curx{j}, 'rl');
        end;
        
        % Now, finally, the optimization itself
%         for j=(L(i)+1):-1:1
        for j=1:(L(i)+1)
            Phi1 = phiaf{i}{j};
            A1 = cura{j};
            Phi2 = phiaf{i}{j+1};
            x1 = curx{j};
            
            y_prev = reshape(cury{j}, ryf{i}(j)*n{i}(j)*ryf{i}(j+1), 1);
            
            dir = 1;
            if (j==(L(i)+1)); dir=0; end;
            local_kickrank = kickrank;
            if (last_sweep); local_kickrank=0; end;
            [u,v,max_dx]=local_proj(Phi1,A1,Phi2, x1, tol2/sqrt(sum(L+1)), y_prev,  max_dx, dir, local_kickrank, verb);
            
            if (j<(L(i)+1))
                cr2 = cury{j+1};
                cr2 = reshape(cr2, ryf{i}(j+1), n{i}(j+1)*ryf{i}(j+2));
                cr2 = (v.')*cr2;
                ryf{i}(j+1) = size(u,2);
                max_frank = max(max_frank, ryf{i}(j+1));
                u = reshape(u, ryf{i}(j), n{i}(j), ryf{i}(j+1));
                cury{j} = u;
                cury{j+1} = reshape(cr2, ryf{i}(j+1), n{i}(j+1), ryf{i}(j+2));
                
                phiaf{i}{j+1} = compute_next_Phi(phiaf{i}{j}, u, cura{j}, x1, 'lr');
            else
                u = u*(v.');
                cury{j} = reshape(u, ryf{i}(j), n{i}(j), ryf{i}(j+1));
            end;
        end;
        if (i<d)
%             % The factor is ready. Now, orthogonalize it l-to-r
%             for j=1:L(i)
%                 cr = curx{j};
%                 cr = reshape(cr, rxf{i}(j)*n{i}(j), rxf{i}(j+1));
%                 [cr,rv]=qr(cr,0);
%                 cr2 = curx{j+1};
%                 cr2 = reshape(cr2, rxf{i}(j+1), n{i}(j+1)*rxf{i}(j+2));
%                 cr2 = rv*cr2;
%                 rxf{i}(j+1) = size(cr, 2);
%                 cr = reshape(cr, rxf{i}(j), n{i}(j), rxf{i}(j+1));
%                 curx{j} = cr;
%                 curx{j+1} = reshape(cr2, rxf{i}(j+1), n{i}(j+1), rxf{i}(j+2));
%                 
%                 phiaf{i}{j+1} = compute_next_Phi(phiaf{i}{j}, cr, cura{j}, cr, 'lr');
%                 phiyf{i}{j+1} = compute_next_Phi(phiyf{i}{j}, cr, [], cury{j}, 'lr');
%             end;
            % First L cells are all what we need for factor
            yf{i}(1:L(i)) = cury(1:L(i));
            % The last one should go to the core. Moreover, we need xc1->xc2 orth
            cr = permute(cury{L(i)+1}, [2, 1, 3]); % now rc1,rt,rc2
            cr = reshape(cr, ryc(i)*ryf{i}(L(i)+1), ryc(i+1));
            [cr,rv]=qr(cr,0);
            cr2 = yf{i+1}{L(i+1)+1};
            cr2 = reshape(cr2, ryc(i+1), ryf{i+1}(L(i+1)+1)*ryc(i+2));
            cr2 = rv*cr2;
            ryc(i+1) = size(cr,2);
            % Update Factor-sizes too
            n{i+1}(L(i+1)+1)=ryc(i+1);
            ryf{i}(L(i)+2)=ryc(i+1);
            cr = reshape(cr, ryc(i), ryf{i}(L(i)+1), ryc(i+1));
            yf{i}{L(i)+1} = cr;
            yf{i+1}{L(i+1)+1} = reshape(cr2, ryc(i+1), ryf{i+1}(L(i+1)+1), ryc(i+2));
            
            % Update acp, ycp
            % core-projected matrix block
            acp{i} = permute(phiaf{i}{L(i)+1}, [2, 1, 3]);
            acp{i} = reshape(acp{i}, raf{i}(L(i)+1), ryf{i}(L(i)+1)*rxf{i}(L(i)+1));
            acp{i} = reshape(permute(af{i}{L(i)+1}, [1,3,2]), rac(i)*rac(i+1), raf{i}(L(i)+1)) * acp{i};
            acp{i} = reshape(acp{i}, rac(i), rac(i+1), ryf{i}(L(i)+1), rxf{i}(L(i)+1));
            acp{i} = permute(acp{i}, [1, 3, 4, 2]);
%             % core-projected rhs block
%             ycp{i} = phiyf{i}{L(i)+1}.';
%             ycp{i} = reshape(permute(yf{i}{L(i)+1}, [1,3,2]), ryc(i)*ryc(i+1), ryf{i}(L(i)+1)) * ycp{i};
%             ycp{i} = reshape(ycp{i}, ryc(i), ryc(i+1), rxf{i}(L(i)+1));
%             ycp{i} = permute(ycp{i}, [1, 3, 2]);
            
            phiac{i+1} = compute_next_Phi(phiac{i}, cr, acp{i}, xf{i}{L(i)+1}, 'lr');
        else
            yf{i}(1:L(i)) = cury(1:L(i));
            yf{i}{L(i)+1} = permute(cury{L(i)+1}, [2, 1, 3]); % core block
        end;
    end;
    
    
    % Residual check, etc
    if (verb>0)
        fprintf('=mvrk2= sweep %d, max_dx: %3.3e, mrank_c: %d, mrank_f: %d\n', swp, max_dx, max(ryc), max_frank);
    end;
    
    if (last_sweep)
        break;
    end;
    
    if (kickrank<0)
        kickrank=kickrank-1;
    end;
    
%     if (max_dx<tol)&&(kickrank<=-als_iters)
%         last_sweep=true;
%         kickrank=0;
%     end;

    if (max_dx<tol)
        kickrank = 0;
        last_sweep=true;
    end;

%     if (max_dx_prev<=tol*als_tol_high)
%         if (max_dx>max_dx_prev); regurg_cnt=regurg_cnt+1; fprintf('---- Regurgitation %d\n', regurg_cnt); end;
%         if ((regurg_cnt>als_iters)||(max_dx<tol)); last_sweep=true; end;
%         if ((regurg_cnt>0)||(max_dx<=tol*als_tol_low))&&(kickrank>=0); kickrank = -1; end;
%     end;
    
    if (swp==nswp-1)
        last_sweep=true;
    end;
    
    max_dx_prev = max_dx;
    max_dx = 0;
    max_frank = 0;
end;

% Stuff back
for i=1:d
    y.tuck{i} = cell2core(y.tuck{i}, yf{i}(1:L(i)));
    yc{i} = yf{i}{L(i)+1};    
end;
y.core = cell2core(y.core, yc);
y.dphys = d;
% x.sz = n;
end


function [u,v,max_dx]=local_proj(Phi1,A1,Phi2, x1, tol, y_prev,  max_dx, dir, kickrank, verb)
ry1 = size(Phi1,1);
n = size(A1,2);
m = size(A1,3);
ry2 = size(Phi2,1);
ra1 = size(Phi1,2);
ra2 = size(Phi2,2);

rx1 = size(x1,1);
rx2 = size(x1,3);
if (dir>0)
    y_new = reshape(Phi1, ry1*ra1, rx1);
    y_new = y_new*reshape(x1, rx1, m*rx2);
    y_new = reshape(y_new, ry1, ra1, m, rx2);
    y_new = permute(y_new, [2, 3, 1, 4]);
    y_new = reshape(y_new, ra1*m, ry1*rx2);
    A1 = permute(A1, [2, 4, 1, 3]);
    A1 = reshape(A1, n*ra2, ra1*m);
    y_new = A1*y_new;
    y_new = reshape(y_new, n, ra2, ry1, rx2);
    y_new = permute(y_new, [3, 1, 2, 4]);
    y_new = reshape(y_new, ry1*n, ra2*rx2);
    if (kickrank>0); y_save = y_new; end;
    y_new = y_new*(reshape(Phi2, ry2, ra2*rx2).');
else
    y_new = reshape(Phi2, ry2*ra2, rx2);
    y_new = reshape(x1, rx1*m, rx2)*(y_new.');
    y_new = reshape(y_new, rx1, m, ry2, ra2);
    y_new = permute(y_new, [2, 4, 1, 3]);
    y_new = reshape(y_new, m*ra2, rx1*ry2);
    A1 = reshape(A1, ra1*n, m*ra2);
    y_new = A1*y_new;
    y_new = reshape(y_new, ra1, n, rx1, ry2);
    y_new = permute(y_new, [1, 3, 2, 4]);
    y_new = reshape(y_new, ra1*rx1, n*ry2);
    if (kickrank>0); y_save = y_new; end;
    y_new = reshape(Phi1, ry1, ra1*rx1)*y_new;
end;


dx = norm(y_new(:)-y_prev(:))/norm(y_new(:));
max_dx = max(max_dx, dx);


% Truncation
if (dir>=0) % left-to-right
    y_new = reshape(y_new, ry1*n, ry2);
else
    y_new = reshape(y_new, ry1, n*ry2);
end;

if (kickrank>=0)
[u,s,v]=svd(y_new, 'econ');
s = diag(s);
    
r = my_chop2(s, tol*norm(s));
    
r = min(r, numel(s));
  
else
    if (dir>=0)
        [u,v]=qr(y_new, 0);
        v = v';
        r = size(u,2);
        s = ones(r,1);
    else
        [v,u]=qr(y_new.', 0);        
        u = u.';
        v = conj(v);
        r = size(v,2);
        s = ones(r,1);        
    end;    
end;

if (verb>1)
    fprintf('=mvrk2= dir %d, dx: %3.3e, r: %d\n', dir, dx, r);
end;
        
    if (dir>0) % left-to-right, kickrank, etc
        u = u(:,1:r);
        v = conj(v(:,1:r))*diag(s(1:r));
        
        if (kickrank>0)
            % Smarter kick: low-rank PCA in residual
            leftresid = reshape(u*v.', ry1*n, ry2);
            leftresid = [leftresid, -y_save];
            
            uk = uchol(leftresid.', kickrank*2);
            uk = uk(:,size(uk,2):-1:max(size(uk,2)-kickrank+1,1));
            
            [u,rv]=qr([u,uk], 0);
            radd = size(uk,2);
            v = [v, zeros(ry2, radd)];
            v = v*(rv.');
        end;
    elseif (dir<0) % right-to-left
        u = u(:,1:r)*diag(s(1:r));
        v = conj(v(:,1:r));
        
        if (kickrank>0)
            % Smarter kick: low-rank PCA in residual
            rightresid = reshape(u*v.', ry1, n*ry2);
            rightresid = [rightresid; -y_save];
            uk = uchol(rightresid, kickrank*2);
            uk = uk(:,size(uk,2):-1:max(size(uk,2)-kickrank+1,1));
            
            [v,rv]=qr([v,uk], 0);
            radd = size(uk,2);
            u = [u, zeros(ry1, radd)];
            u = u*(rv.');
        end;
    else
        % Just stuff back the last core
        u = u(:,1:r);
        v = conj(v(:,1:r))*diag(s(1:r));
    end;

end




function [Phi] = compute_next_Phi(Phi_prev, x, A, y, direction)
% Performs the recurrent Phi (or Psi) matrix computation
% Phi = Phi_prev * (x'Ay).
% If direction is 'lr', computes Psi
% if direction is 'rl', computes Phi
% A can be empty, then only x'y is computed.

if (strcmp(direction, 'rl'))
  % Revert ranks to perform the right-to-left recursion
  x = permute(x, [3, 2, 1]);
  y = permute(y, [3, 2, 1]);
  if (~isempty(A))
    A = permute(A, [4, 2, 3, 1]);
  end
end

rx1 = size(x,1); n = size(x,2); rx2 = size(x,3);
ry1 = size(y,1); m = size(y,2); ry2 = size(y,3);
if (~isempty(A))
  ra1 = size(A,1); ra2 = size(A,4);
else
  ra1 = 1; ra2 = 1;
end

Phi = reshape(Phi_prev, [rx1*ra1, ry1]);
y = reshape(y, [ry1, m*ry2]);
Phi = Phi*y;	% complexity §\mcommentfont$\mathcal{O}(n  r_x r_A r_y^2)$§
Phi = reshape(Phi, [rx1, ra1, m, ry2]);
Phi = permute(Phi, [2, 3, 1, 4]);
if (~isempty(A))
  Phi = reshape(Phi, [ra1*m, rx1*ry2]);
  A = permute(A, [4, 2, 1, 3]);
  A = reshape(A, [ra2*n, ra1*m]);
  Phi = A*Phi;	% complexity §\mcommentfont$\mathcal{O}(n^2  r_x r_A^2 r_y)$§
  Phi = reshape(Phi, [ra2, n, rx1, ry2]);
end
Phi = permute(Phi, [3, 2, 1, 4]);
Phi = reshape(Phi, [rx1*n, ra2*ry2]);
x = reshape(x, [rx1*n, rx2]);
Phi = (x')*Phi;	% complexity §\mcommentfont$\mathcal{O}(n  r_x^2 r_A r_y)$§
if (~isempty(A))
  Phi = reshape(Phi, [rx2, ra2, ry2]);
end
end

