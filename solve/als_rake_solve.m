function [x]=als_rake_solve(a,y,tol,varargin)
%ALS-MINRES method for the solution of linear systems in QTT-Tucker format
%   [X]=ALS_RAKE_SOLVE(A,Y,TOL,VARARGIN) Attempts to solve the linear
%   system A*X = Y with accuracy EPS using the AMR iteration.
%   Matrix A has to be given in the QTT-Tucker, right-hand side Y should be
%   given in the QTT-Tucker format also. Options are provided in form
%   'PropertyName1',PropertyValue1,'PropertyName2',PropertyValue2 and so
%   on. The parameters are set to default (in brackets in the following) 
%   The list of option names and default values are:
%       o x0 - initial approximation [tensor of all ones] 
%       o nswp - maximal number of AMR sweeps [10]
%       o verb - verbosity level, 0-silent, 1-sweep info, 2-block info [1]
%       o kickrank - rank-increasing parameter [5]
%       o max_full_size - maximal size of the local matrix to full solver 
%       [2500]
%       o local_prec: Local preconditioner, '' - none
%         'cjacobi' - block-Jacobi, diagonal over rank indices ['']
%       o local_iters - number of local gmres restarts [2]
%       o local_restart - dimension of local gmres [40]
%       o trunc_norm - truncation and stopping tolerance: 'resid' - the
%       residual norm is used, 'fro' - Frobenius (L2)  ['resid']
%       o resid_damp - gap between the local solver and truncation. Larger
%       values help to filter some noise, but require more local
%       iterations   [2]
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

Au = [];

% Inner parameters
max_full_size=2500;

resid_damp = 2; % Truncation error to true residual treshold
als_tol_high = 10;
als_tol_low = 4;
als_iters=2;

nswp=10;
local_restart=40;
local_iters=2;

local_prec = '';
% local_prec_char = 0;
% local_prec = 'jacobi';

rmax=Inf;
trunc_norm = 'residual';
% trunc_norm_char = 1;
% trunc_norm = 'fro';

% local_solver = 'gmres';
% local_solver = 'pcg';

verb=1;
kickrank = 5;
x=[];

nswp_f = 1;
nswp_c = 1;

for i=1:2:length(varargin)-1
    switch lower(varargin{i})
        case 'nswp'
            nswp=varargin{i+1};
        case 'rmax'
            rmax=varargin{i+1};
        case 'x0'
            x=varargin{i+1};
        case 'verb'
            verb=varargin{i+1};
        case 'local_prec'
            local_prec=varargin{i+1};
        case 'local_restart'
            local_restart=varargin{i+1};
        case 'local_iters'
            local_iters=varargin{i+1};
%         case 'local_solver'
%             local_solver=varargin{i+1};            
        case 'kickrank'
            kickrank=varargin{i+1};
        case 'als_tol_high'
            als_tol_high=varargin{i+1};            
        case 'als_tol_low'
            als_tol_low=varargin{i+1};                        
        case 'als_iters'
            als_iters=varargin{i+1};                                    
        case  'max_full_size'
            max_full_size=varargin{i+1};
%         case 'step_dpow'
%             step_dpow=varargin{i+1};
%         case 'min_dpow'
%             min_dpow=varargin{i+1};
        case 'resid_damp'
            resid_damp = varargin{i+1};
        case 'trunc_norm'
            trunc_norm = varargin{i+1};
%         case 'bot_conv'
%             bot_conv=varargin{i+1};
%         case 'top_conv'
%             top_conv=varargin{i+1};          
            
        otherwise
            error('Unrecognized option: %s\n',varargin{i});
    end
end

tol2 = tol;

d = y.core.d;
yc = core2cell(y.core);
ryc = y.core.r;
n = cell(d,1);
L = zeros(d,1);
yf = cell(d,1);
ryf = cell(d,1);
for i=1:d
    n{i} = y.tuck{i}.n;
    L(i) = y.tuck{i}.d;
    ryf{i} = y.tuck{i}.r;
    ryf{i}(L(i)+2) = ryc(i+1);
    yf{i} = cell(L(i)+1,1);
    yf{i}(1:L(i)) = core2cell(y.tuck{i});
    yf{i}{L(i)+1} = yc{i};
end;

ac = core2cell(a.core);
rac = a.core.r;
af = cell(d,1);
raf = cell(d,1);
for i=1:d
    raf{i} = a.tuck{i}.r;
    raf{i}(L(i)+2) = rac(i+1);    
    af{i} = cell(L(i)+1,1);
    af{i}(1:L(i)) = core2cell(a.tuck{i});
    af{i}{L(i)+1} = ac{i};
end;

xf = cell(d,1);
rxf = cell(d,1);
if (isempty(x))
    x = qtt_tucker;
    xc = tt_ones(1, d);
    x.core = xc;
    rxc = xc.r;
    xc = core2cell(xc);
    x.tuck = cell(d,1);
    for i=1:d
        rxf{i} = ones(L(i)+2,1);
        n{i}(L(i)+1) = 1;
        xf{i} = cell(L(i)+1,1);
        x.tuck{i} = tt_ones(n{i}(1:L(i)));
        xf{i}(1:L(i)) = core2cell(x.tuck{i});
        xf{i}{L(i)+1}=xc{i};
    end;
else
    xc = core2cell(x.core);
    rxc = x.core.r;
    for i=1:d
        rxf{i} = x.tuck{i}.r;
        rxf{i}(L(i)+2)=rxc(i+1);
        n{i}(L(i)+1) = rxc(i);
        xf{i} = cell(L(i)+1,1);
        xf{i}(1:L(i)) = core2cell(x.tuck{i});
        xf{i}{L(i)+1}=xc{i};
    end;    
end;

max_dx = 0;
max_res = 0;

phiaf = cell(d,1);
phiyf = cell(d,1);
for i=1:d
    phiaf{i} = cell(L(i)+2,1);
    phiaf{i}{1}=1;
    phiyf{i} = cell(L(i)+2,1);
    phiyf{i}{1}=1;    
end;
phiac = cell(d+1,1); phiac{1}=1; phiac{d+1}=1;
phiyc = cell(d+1,1); phiyc{1}=1; phiyc{d+1}=1;

acp = cell(d,1);
ycp = cell(d,1);


max_res_prev = Inf;
max_dx_prev = Inf;
max_frank = 0;
regurg_cnt = 0;
last_sweep = false;

for swp=1:nswp
    
    % Orthog + phi, Factors
    for i=1:d
        n{i}(L(i)+1)=rxc(i);
        rxf{i}(L(i)+2)=rxc(i+1);        
        % Permute core blocks: rc1 is mode index now
        xf{i}{L(i)+1} = permute(xf{i}{L(i)+1}, [2,1,3]); % rxt,rxc1,rxc2
        for j=1:L(i)
            cr = xf{i}{j};
            cr = reshape(cr, rxf{i}(j)*n{i}(j), rxf{i}(j+1));
            [cr,rv]=qr(cr,0);
            cr2 = xf{i}{j+1};
            cr2 = reshape(cr2, rxf{i}(j+1), n{i}(j+1)*rxf{i}(j+2));
            cr2 = rv*cr2;
            rxf{i}(j+1) = size(cr, 2);
            cr = reshape(cr, rxf{i}(j), n{i}(j), rxf{i}(j+1));
            xf{i}{j} = cr;
            xf{i}{j+1} = reshape(cr2, rxf{i}(j+1), n{i}(j+1), rxf{i}(j+2));
            
            phiaf{i}{j+1} = compute_next_Phi(phiaf{i}{j}, cr, af{i}{j}, cr, 'lr');
            phiyf{i}{j+1} = compute_next_Phi(phiyf{i}{j}, cr, [], yf{i}{j}, 'lr');
        end;
        % core-projected matrix block
        acp{i} = permute(phiaf{i}{L(i)+1}, [2, 1, 3]);
        acp{i} = reshape(acp{i}, raf{i}(L(i)+1), rxf{i}(L(i)+1)*rxf{i}(L(i)+1));
        acp{i} = reshape(permute(af{i}{L(i)+1}, [1,3,2]), rac(i)*rac(i+1), raf{i}(L(i)+1)) * acp{i};
        acp{i} = reshape(acp{i}, rac(i), rac(i+1), rxf{i}(L(i)+1), rxf{i}(L(i)+1));
        acp{i} = permute(acp{i}, [1, 3, 4, 2]);
        % core-projected rhs block
        ycp{i} = phiyf{i}{L(i)+1}.';
        ycp{i} = reshape(permute(yf{i}{L(i)+1}, [1,3,2]), ryc(i)*ryc(i+1), ryf{i}(L(i)+1)) * ycp{i};
        ycp{i} = reshape(ycp{i}, ryc(i), ryc(i+1), rxf{i}(L(i)+1));
        ycp{i} = permute(ycp{i}, [1, 3, 2]);
        % return xf{i}{L{i}+1} to the core state
        xf{i}{L(i)+1} = permute(xf{i}{L(i)+1}, [2,1,3]); % rxc1,rxt,rxc2
    end;
    
    % ALS over core
    nc = zeros(d,1);
    for i=1:d
        nc(i)=rxf{i}(L(i)+1);
    end;
    for k=1:nswp_c
%     max_res_c = 0;        
    % Orthog
    for i=1:d-1
        cr = xf{i}{L(i)+1};
        cr = reshape(cr, rxc(i)*nc(i), rxc(i+1));
        [cr,rv]=qr(cr,0);
        cr2 = xf{i+1}{L(i+1)+1};
        cr2 = reshape(cr2, rxc(i+1), nc(i+1)*rxc(i+2));
        cr2 = rv*cr2;
        rxc(i+1) = size(cr,2);
        % Update Factor-sizes too
%         n{i+1}(L(i+1)+1)=rxc(i+1);
%         rxf{i}(L(i)+2)=rxc(i+1);
        cr = reshape(cr, rxc(i), nc(i), rxc(i+1));
        xf{i}{L(i)+1} = cr;
        xf{i+1}{L(i+1)+1} = reshape(cr2, rxc(i+1), nc(i+1), rxc(i+2));
        
        phiac{i+1} = compute_next_Phi(phiac{i}, cr, acp{i}, cr, 'lr');
        phiyc{i+1} = compute_next_Phi(phiyc{i}, cr, [], ycp{i}, 'lr');
    end;
    % ALS
    for i=d:-1:1
        Phi1 = phiac{i};
        A1 = acp{i};
        Phi2 = phiac{i+1};
        
        sol_prev = reshape(xf{i}{L(i)+1}, rxc(i)*nc(i)*rxc(i+1), 1);
        
        dir = -1;
        if (i==1); dir=0; end;
        local_kickrank = kickrank;
        if (last_sweep); local_kickrank=0; end;
        if (verb>1); fprintf('\t core %d\n', i); end;
        [u,v,max_res,max_dx,flg]=local_solve(Phi1,A1,Phi2, phiyc{i},ycp{i},phiyc{i+1}, ...
            tol2/sqrt(sum(L+1))/resid_damp, resid_damp, trunc_norm, sol_prev, ...
            local_prec, local_restart, local_iters, max_full_size, max_res, max_dx, ...
            dir, rmax, local_kickrank, verb);
        
        if (flg>0); fprintf('-warn- local_solve did not converge at cb {%d}\n', i); end;
        
        if (i>1)
            cr2 = xf{i-1}{L(i-1)+1};
            cr2 = reshape(cr2, rxc(i-1)*nc(i-1), rxc(i));
            cr2 = cr2*u;
            rxc(i) = size(v,2);
            % Update Factor-sizes
%             n{i}(L(i)+1)=rxc(i);
%             rxf{i-1}(L(i-1)+2)=rxc(i);
            v = reshape(v.', rxc(i), nc(i), rxc(i+1));
            xf{i}{L(i)+1} = v;
            xf{i-1}{L(i-1)+1} = reshape(cr2, rxc(i-1), nc(i-1), rxc(i));
            
            phiac{i} = compute_next_Phi(phiac{i+1}, v, acp{i}, v, 'rl');
            phiyc{i} = compute_next_Phi(phiyc{i+1}, v, [], ycp{i}, 'rl');            
        else
            v = u*(v.');
            xf{i}{L(i)+1} = reshape(v, rxc(i), nc(i), rxc(i+1));
        end;
    end;
    end;
%     max_res = max(max_res, max_res_c);
    
    
    % Optimization over factors
    for i=1:d
        % ALS L->1, the orth. holds
        % Copy the extended factor data to the current array
        curx = xf{i};
        % Last mode size <- largest core rank
        if (rxc(i)>rxc(i+1))
            % Reshape the core block to be the factor one
            curx{L(i)+1} = permute(curx{L(i)+1}, [2, 1, 3]); % rtx, rxc1, rxc2
            
            cura = af{i};
            % last block has to be convolved with phiac-left
            cura{L(i)+1} = permute(phiac{i}, [1,3,2]);
            cura{L(i)+1} = reshape(cura{L(i)+1}, rxc(i)*rxc(i), rac(i));
            cura{L(i)+1} = cura{L(i)+1}*reshape(af{i}{L(i)+1}, rac(i), raf{i}(L(i)+1)*rac(i+1));
            cura{L(i)+1} = reshape(cura{L(i)+1}, rxc(i), rxc(i), raf{i}(L(i)+1), rac(i+1));
            cura{L(i)+1} = permute(cura{L(i)+1}, [3, 1,2, 4]); % tucker rank is now the left
            
            cury = yf{i};
            cury{L(i)+1} = phiyc{i}*reshape(yf{i}{L(i)+1}, ryc(i), ryf{i}(L(i)+1)*ryc(i+1));
            cury{L(i)+1} = reshape(cury{L(i)+1}, rxc(i), ryf{i}(L(i)+1), ryc(i+1));
            cury{L(i)+1} = permute(cury{L(i)+1}, [2, 1, 3]); % tucker rank is now the left
            
            % The L+2-th phi comes from the core
            phiaf{i}{L(i)+2} = phiac{i+1}; % sizes rxc2,ra2,rxc2. rxc2 - our last "rank" index
            phiyf{i}{L(i)+2} = phiyc{i+1};
            
            n{i}(L(i)+1) = rxc(i);
            rxf{i}(L(i)+2) = rxc(i+1);
        else
            % Reshape the core block to be the factor one
            curx{L(i)+1} = permute(curx{L(i)+1}, [2, 3, 1]); %  rtx, rxc2, rxc1
            
            cura = af{i};
            % last block has to be convolved with phiac-left
            cura{L(i)+1} = permute(phiac{i+1}, [2, 1,3]);
            cura{L(i)+1} = reshape(cura{L(i)+1}, rac(i+1), rxc(i+1)*rxc(i+1));
            cura{L(i)+1} = reshape(af{i}{L(i)+1}, rac(i)*raf{i}(L(i)+1), rac(i+1)) * cura{L(i)+1};
            cura{L(i)+1} = reshape(cura{L(i)+1}, rac(i), raf{i}(L(i)+1), rxc(i+1), rxc(i+1));
            cura{L(i)+1} = permute(cura{L(i)+1}, [2, 3,4, 1]); % tucker rank is now the left
            
            cury = yf{i};
            cury{L(i)+1} = reshape(yf{i}{L(i)+1}, ryc(i)*ryf{i}(L(i)+1), ryc(i+1)) * (phiyc{i+1}.');
            cury{L(i)+1} = reshape(cury{L(i)+1}, ryc(i), ryf{i}(L(i)+1), rxc(i+1));
            cury{L(i)+1} = permute(cury{L(i)+1}, [2, 3, 1]); % tucker rank is now the left
            
            % The L+2-th phi comes from the core
            phiaf{i}{L(i)+2} = phiac{i}; % sizes rxc2,ra2,rxc2. rxc2 - our last "rank" index
            phiyf{i}{L(i)+2} = phiyc{i};    
            
            n{i}(L(i)+1) = rxc(i+1);
            rxf{i}(L(i)+2) = rxc(i);
        end;
        
        for k=1:nswp_f
%         max_res_f = 0;
        % First, orthogonality L->1
        for j=(L(i)+1):-1:2
            cr = curx{j};
            cr = reshape(cr, rxf{i}(j), n{i}(j)*rxf{i}(j+1));
            [cr,rv]=qr(cr.',0);
            cr2 = curx{j-1};
            cr2 = reshape(cr2, rxf{i}(j-1)*n{i}(j-1), rxf{i}(j));
            cr2 = cr2*(rv.');
            rxf{i}(j) = size(cr, 2);
            cr = reshape(cr.', rxf{i}(j), n{i}(j), rxf{i}(j+1));
            curx{j} = cr;
            curx{j-1} = reshape(cr2, rxf{i}(j-1), n{i}(j-1), rxf{i}(j));
            
            phiaf{i}{j} = compute_next_Phi(phiaf{i}{j+1}, cr, cura{j}, cr, 'rl');
            phiyf{i}{j} = compute_next_Phi(phiyf{i}{j+1}, cr, [], cury{j}, 'rl');
        end;
        
        % Now, finally, the optimization itself
%         for j=(L(i)+1):-1:1
        for j=1:(L(i)+1)
            Phi1 = phiaf{i}{j};
            A1 = cura{j};
            Phi2 = phiaf{i}{j+1};
            
            sol_prev = reshape(curx{j}, rxf{i}(j)*n{i}(j)*rxf{i}(j+1), 1);
            
            dir = 1;
            if (j==(L(i)+1)); dir=0; end;
            local_kickrank = kickrank;
            if (last_sweep); local_kickrank=0; end;
            if (verb>1); fprintf('\t factor {%d,%d}\n', i, j); end;
            [u,v,max_res,max_dx,flg]=local_solve(Phi1,A1,Phi2, phiyf{i}{j},cury{j},phiyf{i}{j+1}, ...
                tol2/sqrt(sum(L+1))/resid_damp, resid_damp, trunc_norm, sol_prev, ...
                local_prec, local_restart, local_iters, max_full_size, max_res, max_dx, ...
                dir, rmax, local_kickrank, verb);
            
            if (flg>0); fprintf('-warn- local_solve did not converge at fb {%d}{%d}\n', i, j); end;
            
            if (j<(L(i)+1))
                cr2 = curx{j+1};
                cr2 = reshape(cr2, rxf{i}(j+1), n{i}(j+1)*rxf{i}(j+2));
                cr2 = (v.')*cr2;
                rxf{i}(j+1) = size(u,2);
                max_frank = max(max_frank, rxf{i}(j+1));
                u = reshape(u, rxf{i}(j), n{i}(j), rxf{i}(j+1));
                curx{j} = u;
                curx{j+1} = reshape(cr2, rxf{i}(j+1), n{i}(j+1), rxf{i}(j+2));
                
                phiaf{i}{j+1} = compute_next_Phi(phiaf{i}{j}, u, cura{j}, u, 'lr');
                phiyf{i}{j+1} = compute_next_Phi(phiyf{i}{j}, u, [], cury{j}, 'lr');
            else
                u = u*(v.');
                curx{j} = reshape(u, rxf{i}(j), n{i}(j), rxf{i}(j+1));
            end;
        end;
        end;
%         max_res = max(max_res, max_res_f);
        
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
            xf{i}(1:L(i)) = curx(1:L(i));
            % The last one should go to the core. Moreover, we need xc1->xc2 orth
            if (rxc(i)>rxc(i+1))
                cr = permute(curx{L(i)+1}, [2, 1, 3]); % now rc1,rt,rc2
            else 
                cr = permute(curx{L(i)+1}, [3, 1, 2]); % now rc1,rt,rc2
            end;
            cr = reshape(cr, rxc(i)*rxf{i}(L(i)+1), rxc(i+1));
            [cr,rv]=qr(cr,0);
            cr2 = xf{i+1}{L(i+1)+1};
            cr2 = reshape(cr2, rxc(i+1), rxf{i+1}(L(i+1)+1)*rxc(i+2));
            cr2 = rv*cr2;
            rxc(i+1) = size(cr,2);
            % Update Factor-sizes too
%             n{i+1}(L(i+1)+1)=rxc(i+1);
%             rxf{i}(L(i)+2)=rxc(i+1);
            cr = reshape(cr, rxc(i), rxf{i}(L(i)+1), rxc(i+1));
            xf{i}{L(i)+1} = cr;
            xf{i+1}{L(i+1)+1} = reshape(cr2, rxc(i+1), rxf{i+1}(L(i+1)+1), rxc(i+2));
            
            % Update acp, ycp
            % core-projected matrix block
            acp{i} = permute(phiaf{i}{L(i)+1}, [2, 1, 3]);
            acp{i} = reshape(acp{i}, raf{i}(L(i)+1), rxf{i}(L(i)+1)*rxf{i}(L(i)+1));
            acp{i} = reshape(permute(af{i}{L(i)+1}, [1,3,2]), rac(i)*rac(i+1), raf{i}(L(i)+1)) * acp{i};
            acp{i} = reshape(acp{i}, rac(i), rac(i+1), rxf{i}(L(i)+1), rxf{i}(L(i)+1));
            acp{i} = permute(acp{i}, [1, 3, 4, 2]);
            % core-projected rhs block
            ycp{i} = phiyf{i}{L(i)+1}.';
            ycp{i} = reshape(permute(yf{i}{L(i)+1}, [1,3,2]), ryc(i)*ryc(i+1), ryf{i}(L(i)+1)) * ycp{i};
            ycp{i} = reshape(ycp{i}, ryc(i), ryc(i+1), rxf{i}(L(i)+1));
            ycp{i} = permute(ycp{i}, [1, 3, 2]);
            
            phiac{i+1} = compute_next_Phi(phiac{i}, cr, acp{i}, cr, 'lr');
            phiyc{i+1} = compute_next_Phi(phiyc{i}, cr, [], ycp{i}, 'lr');
        else
            xf{i}(1:L(i)) = curx(1:L(i));
            if (rxc(i)>rxc(i+1))
                xf{i}{L(i)+1} = permute(curx{L(i)+1}, [2, 1, 3]); % core block
            else
                xf{i}{L(i)+1} = permute(curx{L(i)+1}, [3, 1, 2]); % core block
            end;
        end;
    end;
    
    
    % Residual check, etc
    
%     x_old = x;
%     for i=1:d
%         x.tuck{i} = cell2core(x.tuck{i}, xf{i}(1:L(i)));
%         xc{i} = xf{i}{L(i)+1};
%     end;
%     x.core = cell2core(x.core, xc);
    
    if (verb>0)
        fprintf('=als_rake_solve= sweep %d, max_dx: %3.3e, max_res: %3.3e, mrank_c: %d, mrank_f: %d\n', swp, max_dx, max_res, max(rxc), max_frank);
    end;
    
    if (last_sweep)
        break;
    end;    
    
    if (kickrank<0)
        kickrank=kickrank-1;
    end;    
    
    if (strcmp(trunc_norm, 'fro'))        
        if (max_dx<tol)&&(kickrank<=-als_iters)
            last_sweep=true;
%             tol2 = tol;
        end;
        if (max_dx_prev<=tol*als_tol_high)
            if (max_dx>max_dx_prev); regurg_cnt=regurg_cnt+1; fprintf('---- Regurgitation %d\n', regurg_cnt); end;           
%             kickrank = -1;
%             tol2=tol/4;
        end;
        if ((regurg_cnt>0)||(max_dx<=tol*als_tol_low))&&(kickrank>0); kickrank=-1; end;
    else
%         for i=1:d
%             x.tuck{i} = cell2core(x.tuck{i}, xf{i}(1:L(i)));
%             xc{i} = xf{i}{L(i)+1};
%         end;
%         x.core = cell2core(x.core, xc);
%         x.dphys = d;
%         Au = mvrk2(a, x, tol, 'y0', Au, 'verb', 0);
%         real_res = norm(Au-y)/norm(y);
%         fprintf('=als_rake_solve= sweep %d, \t\t\t real_res: %3.3e\n', swp, real_res);
        if (max_res<tol)
            break;
%             kickrank = 0;
%             last_sweep = true;
        end;
%         if (max_res<tol)&&(kickrank<=-als_iters)
%             last_sweep=true;
%         end;
%         if (max_res<=tol*als_tol_high)
%             if (max_res>max_res_prev); regurg_cnt=regurg_cnt+1; fprintf('---- Regurgitation %d\n', regurg_cnt); end;
%             if ((regurg_cnt>als_iters)||(max_res<tol)); last_sweep=true; end;
% %             if ((regurg_cnt>0)||(max_res<=tol*als_tol_low))&&(kickrank>0); kickrank=-1; end;
%         end;        
    end;
    
    if (swp==nswp-1)||(regurg_cnt>2)
        last_sweep=true;
%         tol2 = tol;
    end;
    
    max_res_prev = max_res;
    max_dx_prev = max_dx;
    
    max_res = 0;
    max_dx = 0;
    max_frank=0;
end;

% Stuff back
for i=1:d
    x.tuck{i} = cell2core(x.tuck{i}, xf{i}(1:L(i)));
    xc{i} = xf{i}{L(i)+1};    
end;
x.core = cell2core(x.core, xc);
x.dphys = d;
% x.sz = n;
end



function [y]=bfun3(Phi1,B1,Phi2, x)
% Computes (Phi1 * B1 * Phi2)*x
% Phi1 is of sizes ry1, rB1, rx1
% B1 is of sizes rB1, k1, m1, rB2
% Phi2 is of sizes ry2, rB2, rx2
ry1 = size(Phi1,1); ry2 = size(Phi2,1);
rx1 = size(Phi1,3); rx2 = size(Phi2,3);
rb1=size(B1,1); rb2=size(B1,4); 
m1 = size(B1,3);
k1 = size(B1,2);

y = reshape(x, rx1, m1*rx2);
Phi1 = reshape(Phi1, ry1*rb1, rx1);
y = Phi1*y; % size ry1*rb1,m1*rx2 % cplx rb*rx^3*m^2
y = reshape(y, ry1, rb1*m1, rx2);
y = permute(y, [2, 1, 3]);
y = reshape(y, rb1*m1, ry1*rx2);
B1 = permute(B1, [2, 4, 1, 3]);
B1 = reshape(B1, k1*rb2, rb1*m1);
y = B1*y; % size k1*rb2, ry1*rx2 % cplx rb^2*rx^2*n^3
y = reshape(y, k1, rb2, ry1, rx2);
y = permute(y, [2, 4, 3, 1]);
y = reshape(y, rb2*rx2, ry1*k1);
Phi2 = reshape(Phi2, ry2, rb2*rx2);
y = Phi2*y; % size ry2, ry1*k1 % cplx rb*rx^3*n^2
y = y.';
y = reshape(y, ry1*k1*ry2, 1);
end


function [u,v,max_res,max_dx,flg]=local_solve(Phi1,A1,Phi2, phiy1,y1,phiy2, tol, resid_damp, trunc_norm, sol_prev, local_prec, local_restart, local_iters, max_full_size, max_res, max_dx, dir, rmax, kickrank, verb)
rx1 = size(Phi1,1);
n = size(A1,2);
rx2 = size(Phi2,1);
ra1 = size(Phi1,2);
ra2 = size(Phi2,2);

ry1 = size(phiy1,2);
ry2 = size(phiy2,2);
if (dir>0)
    rhs = reshape(y1, ry1, n*ry2);
    rhs = phiy1*rhs;
    rhs = reshape(rhs, rx1*n, ry2);
    if (kickrank>0); y_save = rhs; end;
    rhs = rhs*(phiy2.');
else
    rhs = phiy2.';
    rhs = reshape(y1, ry1*n, ry2)*rhs;
    rhs = reshape(rhs, ry1, n*rx2);
    if (kickrank>0); y_save = rhs; end;
    rhs = phiy1*rhs;
end;
rhs = rhs(:);

norm_rhs = norm(rhs, 'fro');

if (rx1*n*rx2<max_full_size)
    %      |     |    |
    % B = Phi1 - A1 - Phi2
    %      |     |    |
    B = reshape(permute(Phi1, [1, 3, 2]), rx1*rx1, ra1);
    B = B*reshape(A1, ra1, n*n*ra2);
    B = reshape(B, rx1, rx1, n, n, ra2);
    B = permute(B, [1, 3, 2, 4, 5]);
    B = reshape(B, rx1*n*rx1*n, ra2);
    B = B*reshape(permute(Phi2, [2, 1, 3]), ra2, rx2*rx2);
    B = reshape(B, rx1*n, rx1*n, rx2, rx2);
    B = permute(B, [1, 3, 2, 4]);
    B = reshape(B, rx1*n*rx2, rx1*n*rx2);
    
    if (norm_rhs==0.0)
        % Ground state
        res_prev = norm(B*sol_prev);
        if (res_prev>tol)
            B2 = B+eye(rx1*n*rx2);
            sol_prev2 = sol_prev;
            for it=1:local_restart
                sol = B2 \ sol_prev2;
                sol = sol/norm(sol);
                res_new = norm(B*sol);
                if (strcmp(trunc_norm, 'fro'))
                    if (norm(sol-sol_prev2)<tol); break; end;
                else
                    if (res_new<tol); break; end;
                end;
                sol_prev2 = sol;
            end;
            flg = 0;
        else
            sol = sol_prev;
            res_new = res_prev;
            flg=0;
        end;
    else
        res_prev = norm(B*sol_prev-rhs)/norm_rhs;
        if (res_prev>tol)
            sol = B \ rhs;
            % If the system was ill-conditioned
            %         [sol,flg] = gmres(B, rhs, local_restart, real_tol, 2, [], [], sol);
            res_new = norm(B*sol-rhs)/norm_rhs;
            flg=0;
            if (res_new>tol); flg=1; end;
        else
            sol = sol_prev;
            res_new = res_prev;
            flg=0;
        end;
    end;
else
    if (norm_rhs==0.0)
        res_prev = norm(bfun3(Phi1, A1, Phi2, sol_prev));
        if (res_prev>tol)
            trunc_norm_char = 1;
            local_prec_char = 0;
            if ((strcmp(local_prec, 'cjacobi'))); local_prec_char = 1; end;
            if (strcmp(local_prec, 'ljacobi')); local_prec_char = 2;  end;
            if (strcmp(local_prec, 'rjacobi')); local_prec_char = 3;  end;
            
            Psi1 = zeros(rx1, rx1, ra1+1);
            Psi1(:,:,1:ra1) = permute(Phi1, [1,3,2]);
            Psi1(:,:,ra1+1) = eye(rx1);
            B1 = zeros(n, n, ra1+1, ra2+1);
            B1(:,:,1:ra1,1:ra2) = permute(A1, [2,3,1,4]);
            B1(:,:,ra1+1,ra2+1) = eye(n);
            B1 = permute(B1, [3,1,2,4]);
            Psi2 = zeros(rx2,rx2,ra2+1);
            Psi2(:,:,1:ra2) = permute(Phi2, [1,3,2]);
            Psi2(:,:,ra2+1) = eye(rx2);
            Psi2 = permute(Psi2, [2,3,1]);
            
            sol_prev2 = sol_prev;
            for it=1:local_restart
                sol = solve3d_2(Psi1, B1, Psi2, sol_prev2, tol, trunc_norm_char, sol_prev2, local_prec_char, local_restart, local_iters, 1);
                sol = sol/norm(sol);
                res_new = norm(bfun3(Phi1, A1, Phi2, sol));
                if (strcmp(trunc_norm, 'fro'))
                    if (norm(sol-sol_prev2)<tol); break; end;
                else
                    if (res_new<tol); break; end;
                end;
                sol_prev2 = sol;
            end;
            
            flg=0;
            if (res_new>tol); flg=1; end;
        else
            sol = sol_prev;
            res_new = res_prev;
            flg = 0;
        end;
    else
        res_prev = norm(bfun3(Phi1, A1, Phi2, sol_prev) - rhs)/norm_rhs;
        if (res_prev>tol)
            trunc_norm_char = 1;
%             if (strcmp(trunc_norm, 'fro')); trunc_norm_char = 0; end;
            local_prec_char = 0;
            if ((strcmp(local_prec, 'cjacobi'))); local_prec_char = 1; end;
            if (strcmp(local_prec, 'ljacobi')); local_prec_char = 2;  end;
            if (strcmp(local_prec, 'rjacobi')); local_prec_char = 3;  end;
            if (res_prev>1)
                sol_prev = zeros(rx1*n*rx2, 1);
                res_prev = 1;
            end;
            
%             sol = solve3d(permute(Phi1,[1,3,2]), A1, permute(Phi2, [1,3,2]), rhs, tol, trunc_norm_char, sol_prev, local_prec_char, local_restart, local_iters, 1);
            sol = solve3d_2(permute(Phi1,[1,3,2]), A1, permute(Phi2, [3,2,1]), rhs, tol, trunc_norm_char, sol_prev, local_prec_char, local_restart, local_iters, 0);
           
            res_new = norm(bfun3(Phi1, A1, Phi2, sol) - rhs)/norm_rhs;
            flg=0;
            if (res_new>tol); flg=1; end;
        else
            sol = sol_prev;
            res_new = res_prev;
            flg = 0;
        end;
    end;
end;

if (res_prev/res_new<resid_damp)&&(res_new>tol)
    fprintf('--warn-- the residual damp was smaller than in the truncation\n');
end;

% if (flg>0)&&(kickrank==0)
%     sol = sol_prev;
%     res_new = res_prev;
% end;

dx = norm(sol-sol_prev)/norm(sol);
max_dx = max(max_dx, dx);

max_res = max(max_res, res_prev);


if (norm_rhs==0.0)
    norm_rhs=1;
end;

% Truncation
if (dir>=0) % left-to-right
    sol = reshape(sol, rx1*n, rx2);
else
    sol = reshape(sol, rx1, n*rx2);
end;

if (kickrank>=0)
[u,s,v]=svd(sol, 'econ');
s = diag(s);
    
if (strcmp(trunc_norm, 'fro')) % We are happy with L2 truncation (when? but let it be)    
    r = my_chop2(s, tol*resid_damp*norm(s));
else
    % Residual trunc; First, bin-search
    r1 = 1; r2 = numel(s); r = round((r1+r2)/2);
    while (r2-r1>1)
        cursol = u(:,1:r)*diag(s(1:r))*(v(:,1:r)');
        if (rx1*n*rx2<max_full_size)
            res = norm(B*cursol(:)-rhs)/norm_rhs;
        else
            res = norm(bfun3(Phi1, A1, Phi2, cursol)-rhs)/norm_rhs;
        end;
        if (res<max(tol, res_new)*resid_damp)
            r2 = r;
        else
            r1 = r;
        end;
        r = round((r1+r2)/2);
    end;
    r = max(r-1,1);
    % More accurate Linear search
    while (r<=numel(s))
        cursol = u(:,1:r)*diag(s(1:r))*(v(:,1:r)');
        if (rx1*n*rx2<max_full_size)
            res = norm(B*cursol(:)-rhs)/norm_rhs;
        else
            res = norm(bfun3(Phi1, A1, Phi2, cursol)-rhs)/norm_rhs;
        end;
        if (res<max(tol, res_new)*resid_damp)
            break;
        end;
        r = r+1;
    end;
end;
    
r = min(r, numel(s));
r = min(r, rmax);

else
    if (dir>=0)
        [u,v]=qr(sol, 0);
        v = v';
        r = size(u,2);
        s = ones(r,1);
    else
        [v,u]=qr(sol.', 0);        
        u = u.';
        v = conj(v);
        r = size(v,2);
        s = ones(r,1);        
    end;
end;

if (verb>1)
    fprintf('=als_rake_solve= dir %d, dx: %3.3e, res_prev: %3.3e, res_new: %3.3e r: %d\n', dir, dx, res_prev, res_new, r);
end;
        
    if (dir>0) % left-to-right, kickrank, etc
        u = u(:,1:r);
        v = conj(v(:,1:r))*diag(s(1:r));
        
        if (kickrank>0)
            % Smarter kick: low-rank PCA in residual
            % Matrix: Phi1-A{i}, rhs: Phi1-y{i}, sizes rx(i)*n - ra(i+1)            
            leftresid = reshape(Phi1, rx1*ra1, rx1)*reshape(u*v.', rx1, n*rx2);
            leftresid = reshape(leftresid, rx1, ra1*n, rx2);
            leftresid = reshape(permute(leftresid, [2, 1, 3]), ra1*n, rx1*rx2);
            leftresid = reshape(permute(A1, [2,4,1,3]), n*ra2, ra1*n)*leftresid;
            leftresid = reshape(leftresid, n, ra2, rx1, rx2);
            leftresid = reshape(permute(leftresid, [3,1,2,4]), rx1*n, ra2*rx2);
%             leftA = permute(Phi1, [1, 3, 2]);
%             leftA = reshape(leftA, rx1*rx1, ra1);
%             leftA = leftA*reshape(A1, ra1, n*n*ra2);
%             leftA = reshape(leftA, rx1, rx1, n, n, ra2);
%             leftA = permute(leftA, [1, 3, 5, 2, 4]);
%             leftA = reshape(leftA, rx1*n*ra2, rx1*n);
%             leftresid = leftA*reshape(u*v.', rx1*n, rx2);
%             leftresid = reshape(leftresid, rx1*n, ra2*rx2);
            leftresid = [leftresid, -y_save];
%             
%             uk = zeros(rx1*n, min(kickrank,rx1*n));
%             for i=1:min(kickrank, rx1*n)
%                 uk2 = uchol(leftresid.', 2);
%                 uk(:,i) = uk2(:,end);
%                 [uk(:,1:i),rv]=qr(uk(:,1:i), 0);
%                 leftresid = leftA*uk(:,i);
%                 leftresid = reshape(leftresid, rx1*n, ra2);
%             end;
            
% %             [uk,~,~]=svd(leftresid, 'econ');
% %             uk = uk(:,1:min(kickrank, size(uk,2)));
            uk = uchol(leftresid.', kickrank+1);
            uk = uk(:,size(uk,2):-1:max(size(uk,2)-kickrank+1,1));
%             leftresid = leftA*uk;
%             leftresid = reshape(leftresid, rx1*n, ra2*size(uk,2));
%             uk(:,size(uk,2)+1) = uchol(leftresid.', 1);
% %             uk = randn(rx1*n, kickrank);

            [u,rv]=qr([u,uk], 0);
            radd = size(uk,2);
            v = [v, zeros(rx2, radd)];
            v = v*(rv.');
        end;
    elseif (dir<0) % right-to-left
        u = u(:,1:r)*diag(s(1:r));
        v = conj(v(:,1:r));
        
        if (kickrank>0)
            % Smarter kick: low-rank PCA in residual
            % Matrix: Phi1-A{i}, rhs: Phi1-y{i}, sizes rx(i)*n - ra(i+1)
            rightresid = reshape(Phi2, rx2*ra2, rx2)*(reshape(u*v.', rx1*n, rx2).');
            rightresid = reshape(rightresid, rx2, ra2, rx1, n);
            rightresid = reshape(permute(rightresid, [4, 2, 3, 1]), n*ra2, rx1*rx2);
            rightresid = reshape(A1, ra1*n, n*ra2)*rightresid;
            rightresid = reshape(rightresid, ra1, n, rx1, rx2);
            rightresid = reshape(permute(rightresid, [2,4,1,3]), n*rx2, ra1*rx1);            
%             rightA = permute(Phi2, [2, 1, 3]);
%             rightA = reshape(rightA, ra2, rx2*rx2);
%             rightA = reshape(A1, ra1*n*n, ra2)*rightA;
%             rightA = reshape(rightA, ra1, n, n, rx2, rx2);
%             rightA = permute(rightA, [2, 4, 1, 3, 5]);
%             rightA = reshape(rightA, n*rx2*ra1, n*rx2);
%             rightresid = rightA*(reshape(u*v.', rx1, n*rx2).');
%             rightresid = reshape(rightresid, n*rx2, ra1*rx1);
            rightresid = [rightresid, -(y_save.')];
            
%             uk = zeros(n*rx2, min(kickrank,n*rx2));
%             for i=1:min(kickrank, n*rx2)
%                 uk2 = uchol(rightresid.', 2);
%                 uk(:,i) = uk2(:,end);
%                 [uk(:,1:i),rv]=qr(uk(:,1:i), 0);
%                 rightresid = rightA*uk(:,i);
%                 rightresid = reshape(rightresid, n*rx2, ra1);
%             end;

% %             [uk,~,~]=svd(rightresid, 'econ');
% %             uk = uk(:,1:min(kickrank, size(uk,2)));

% %             uk = randn(n*rx2, kickrank);
            uk = uchol(rightresid.', kickrank+1);
            uk = uk(:,size(uk,2):-1:max(size(uk,2)-kickrank+1,1));
            
            [v,rv]=qr([v,uk], 0);
            radd = size(uk,2);
            u = [u, zeros(rx1, radd)];
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

