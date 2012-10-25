function [y]=mvrk(a,x,tol,varargin)
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


% rmax=1000;
verb=1;
kickrank = 2;
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
        y.tuck{i} = tt_rand(n{i}(1:L(i)), L(i), 1);
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

phiaf = cell(d,1);
for i=1:d
    phiaf{i} = cell(L(i)+2,1);
    phiaf{i}{1}=1;
end;
phiac = cell(d+1,1); phiac{1}=1; phiac{d+1}=1;

acp = cell(d,1);
% xcp = cell(d,1);

last_sweep = false;

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
        cr = reshape(cr, ryc(i), nc(i), ryc(i+1));
        yf{i}{L(i)+1} = cr;
        yf{i+1}{L(i+1)+1} = reshape(cr2, ryc(i+1), nc(i+1), ryc(i+2));
        
        phiac{i+1} = compute_next_Phi(phiac{i}, cr, acp{i}, xf{i}{L(i)+1}, 'lr');
    end;
    % ALS
    for i=d-1:-1:1
        Phi1 = phiac{i};
        A1 = acp{i};
        A2 = acp{i+1};
        x1 = xf{i}{L(i)+1};
        x2 = xf{i+1}{L(i+1)+1};
        Phi2 = phiac{i+2};
        
        y_prev = reshape(yf{i}{L(i)+1}, ryc(i)*nc(i),  ryc(i+1));
        y_prev = y_prev*reshape(yf{i+1}{L(i+1)+1}, ryc(i+1), nc(i+1)*ryc(i+2));
        
        dir = -1;
        local_kickrank = kickrank;
        if (last_sweep); local_kickrank=0; end;
        [u,v,max_dx]=local_proj2(Phi1,A1,A2,Phi2, x1, x2, tol2/sqrt(sum(L)), y_prev,  max_dx, dir, local_kickrank, verb);
        
        ryc(i+1) = size(v,2);
        u = reshape(u, ryc(i), nc(i), ryc(i+1));
        v = reshape(v.', ryc(i+1), nc(i+1), ryc(i+2));
        yf{i}{L(i)+1} = u;
        yf{i+1}{L(i+1)+1} = v;
        
        phiac{i+1} = compute_next_Phi(phiac{i+2}, v, A2, x2, 'rl');
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
        
        n{i}(L(i)+1) = ryc(i);
        ryf{i}(L(i)+2) = ryc(i+1);
        
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
        for j=1:(L(i))
            Phi1 = phiaf{i}{j};
            A1 = cura{j};
            A2 = cura{j+1};
            Phi2 = phiaf{i}{j+2};
            x1 = curx{j};
            x2 = curx{j+1};
            
            y_prev = reshape(cury{j}, ryf{i}(j)*n{i}(j), ryf{i}(j+1));
            y_prev = y_prev*reshape(cury{j+1}, ryf{i}(j+1), n{i}(j+1)*ryf{i}(j+2));
            
            dir = 1;
            local_kickrank = kickrank;
            if (last_sweep); local_kickrank=0; end;
            [u,v,max_dx]=local_proj2(Phi1,A1,A2,Phi2, x1,x2, tol2/sqrt(sum(L)), y_prev,  max_dx, dir, local_kickrank, verb);
            
            ryf{i}(j+1) = size(u,2);
            u = reshape(u, ryf{i}(j), n{i}(j), ryf{i}(j+1));
            v = reshape(v.', ryf{i}(j+1), n{i}(j+1), ryf{i}(j+2));
            cury{j} = u;
            cury{j+1} = v;
            
            phiaf{i}{j+1} = compute_next_Phi(phiaf{i}{j}, u, A1, x1, 'lr');
        end;
        if (i<d)
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
        fprintf('=mvrk= sweep %d, max_dx: %3.3e, mrank_c: %d, mrank_f: %d\n', swp, max_dx, max(ryc), max(cell2mat(ryf)));
    end;
    
    if (last_sweep)
        break;
    end;
    
    if (max_dx<tol)
        last_sweep=true;
    end;
    
    if (swp==nswp-1)
        last_sweep=true;
    end;
    
    max_dx = 0;
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


function [u,v,max_dx]=local_proj2(Phi1,A1,A2,Phi2, x1,x2, tol, y_prev,  max_dx, dir, kickrank, verb)
ry1 = size(Phi1,1);
n1 = size(A1,2);
m1 = size(A1,3);
n2 = size(A2,2);
m2 = size(A2,3);
ry3 = size(Phi2,1);
ra1 = size(Phi1,2);
ra2 = size(A2,1);
ra3 = size(Phi2,2);

rx1 = size(x1,1);
rx2 = size(x1,3);
rx3 = size(x2,3);

y_new = reshape(Phi1, ry1*ra1, rx1);
y_new = y_new*reshape(x1, rx1, m1*rx2);
y_new = reshape(y_new, ry1, ra1, m1, rx2);
y_new = permute(y_new, [2, 3, 1, 4]);
y_new = reshape(y_new, ra1*m1, ry1*rx2);
A1 = permute(A1, [2, 4, 1, 3]);
A1 = reshape(A1, n1*ra2, ra1*m1);
y_new = A1*y_new;
y_new = reshape(y_new, n1*ra2*ry1, rx2);
y_new = y_new*reshape(x2, rx2, m2*rx3);
y_new = reshape(y_new, n1, ra2, ry1, m2, rx3);
y_new = permute(y_new, [2, 4, 3, 1, 5]);
y_new = reshape(y_new, ra2*m2, ry1*n1*rx3);
A2 = permute(A2, [2,4,1,3]);
A2 = reshape(A2, n2*ra3, ra2*m2);
y_new = A2*y_new;
y_new = reshape(y_new, n2,ra3,ry1*n1,rx3);
y_new = permute(y_new, [3, 1, 2, 4]);
y_new = reshape(y_new, ry1*n1*n2, ra3*rx3);
y_new = y_new*(reshape(Phi2, ry3, ra3*rx3).');


dx = norm(y_new(:)-y_prev(:))/norm(y_new(:));
max_dx = max(max_dx, dx);


y_new = reshape(y_new, ry1*n1, n2*ry3);

[u,s,v]=svd(y_new, 'econ');
s = diag(s);
    
r = my_chop2(s, tol*norm(s));
    
r = min(r, numel(s));
  
if (verb>1)
    fprintf('=mvrk= dir %d, dx: %3.3e, r: %d\n', dir, dx, r);
end;
        
    if (dir>0) % left-to-right, kickrank, etc
        u = u(:,1:r);
        v = conj(v(:,1:r))*diag(s(1:r));
        
        if (kickrank>0)
            uk = randn(ry1*n1, kickrank);
            
            [u,rv]=qr([u,uk], 0);
            radd = size(uk,2);
            v = [v, zeros(n2*ry3, radd)];
            v = v*(rv.');
        end;
    elseif (dir<0) % right-to-left
        u = u(:,1:r)*diag(s(1:r));
        v = conj(v(:,1:r));
        
        if (kickrank>0)
            uk = randn(n2*ry3, kickrank);
            
            [v,rv]=qr([v,uk], 0);
            radd = size(uk,2);
            u = [u, zeros(ry1*n1, radd)];
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

