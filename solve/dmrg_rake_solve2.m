function [x]=dmrg_rake_solve2(A, y, tol, varargin)
%DMRG-type method for the solution of linear systems in QTT-Tucker format
%   [X]=DMRG_RAKE_SOLVE2(A,Y,TOL,VARARGIN) Attempts to solve the linear
%   system A*X = Y with accuracy EPS using the two-sided DMRG iteration.
%   Matrix A has to be given in the QTT-Tucker, right-hand side Y should be
%   given in the QTT-Tucker format also. Options are provided in form
%   'PropertyName1',PropertyValue1,'PropertyName2',PropertyValue2 and so
%   on. The parameters are set to default (in brackets in the following) 
%   The list of option names and default values are:
%       o x0 - initial approximation [random rank-2 tensor] 
%       o nswp - maximal number of DMRG sweeps [10]
%       o rmax - maximal TT-rank of the solution [1000]
%       o verb - verbosity level, 0-silent, 1-sweep info, 2-block info [1]
%       o kick_rank - stabilization parameter [2]
%       o max_full_size - maximal size of the local matrix to full solver 
%       [2500]
%       o local_prec: Local preconditioner, 'als' - ALS-Richardson
%       iteration, 'selfprec' (Saad selfpreconditioner) ['als']
%       o gmres_iters - number of local gmres restarts [2]
%       o nrestart - dimension of local gmres [25]
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

nswp = 20;
local_format = 'full';
% local_format = 'tt';
max_full_size = 2500;
max_full_size2 = Inf;
nrestart = 40;
gmres_iters = 2;
verb = 1;
kickrank = 2;
% checkrank = 1;
resid_damp_loc = 2;
rmax = Inf;
trunc_norm = 'matrix';

tol2 = tol;

x = [];

for i=1:2:length(varargin)-1
    switch lower(varargin{i})
        case 'nswp'
            nswp=varargin{i+1};
        case 'rmax'
           rmax=lower(varargin{i+1});
        case 'x0'
            x=varargin{i+1};
        case 'verb'
            verb=varargin{i+1};
%         case 'local_prec'
%             local_prec=varargin{i+1};
        case 'nrestart'
            nrestart=varargin{i+1};
        case 'gmres_iters'
            gmres_iters=varargin{i+1};
        case 'kickrank'
            kickrank=varargin{i+1};
        case  'max_full_size'
            max_full_size=varargin{i+1};
        case  'resid_damp'
            resid_damp_loc=varargin{i+1};       
        case  'trunc_norm'
            trunc_norm=varargin{i+1};            
%         case 'prec_compr'
%             prec_compr=varargin{i+1};
%         case 'prec_tol'
%             prec_tol=varargin{i+1};
%         case 'prec_iters'
%             prec_iters=varargin{i+1};
%         case 'use_self_prec'
%             use_self_prec=varargin{i+1};
%         case 'ddpow'
%             ddpow=varargin{i+1};
%         case 'ddrank'
%             ddrank=varargin{i+1};
%         case 'd_pow_check'
%             d_pow_check=varargin{i+1};
%         case 'bot_conv'
%             bot_conv=varargin{i+1};
%         case 'top_conv'
%             top_conv=varargin{i+1};
%         case 'min_dpow'
%             min_dpow=varargin{i+1};
%         case 'min_drank'
%             min_drank=varargin{i+1};
        otherwise
            error('Unrecognized option: %s\n',varargin{i});
    end
end

d = y.dphys; % Physical dim.
yc = y.core;
yf = y.tuck;
Af = A.tuck;
Ac = A.core;

L = zeros(1,d); % Quantics dims
n = zeros(max(L), d); % Physical mode sizes
for i=1:d
    L(i) = yf{i}.d;
    n(1:L(i), i) = yf{i}.n;
end;

if (isempty(x))
    xc = tt_rand(2,d,2);
    xf = cell(d,1);
    for i=1:d
        xf{i} = tt_rand(n(1:L(i),i), L(i), [1;2*ones(L(i),1)]);
    end;
else
    xc = x.core;
    xf = x.tuck;
end;


% Extract ranks; Note that rf(L(i)+1,i) = r_tuck(i)
rcy = yc.r;
rfy = zeros(max(L)+1, d);
for i=1:d
    rfy(1:L(i)+1, i) = yf{i}.r;
end;
rcA = Ac.r;
rfA = zeros(max(L)+1, d);
for i=1:d
    rfA(1:L(i)+1, i) = Af{i}.r;
end;
rcx = xc.r;
rfx = zeros(max(L)+1, d);
for i=1:d
    rfx(1:L(i)+1, i) = xf{i}.r;
end;

% Init phis. Thousands of them... (c)
phcA = cell(d+1,1); phcA{1} = 1; phcA{d+1}=1; % Core, matrix
phcy = cell(d+1,1); phcy{1} = 1; phcy{d+1}=1; % Core, rhs
phfA = cell(d,1); % factors, matrix
phfy = cell(d,1); % factors, rhs
phAfc = cell(d,1); % factor-core, matrix
phyfc = cell(d,1); % factor-core, rhs
for i=1:d
    phfA{i} = cell(L(i)+1,1);
    phfA{i}{1} = 1; phfA{i}{L(i)+1} = 1;
    phfy{i} = cell(L(i)+1,1);
    phfy{i}{1} = 1; phfy{i}{L(i)+1} = 1;
end;

% For random check
% cphcA = cell(d+1,1); cphcA{1} = 1; cphcA{d+1}=1; % Core, matrix
% cphcy = cell(d+1,1); cphcy{1} = 1; cphcy{d+1}=1; % Core, rhs
% cphfA = cell(d,1); % factors, matrix
% cphfy = cell(d,1); % factors, rhs
% cphAfc = cell(d,1); % factor-core, matrix
% cphyfc = cell(d,1); % factor-core, rhs
% for i=1:d
%     cphfA{i} = cell(L(i)+1,1);
%     cphfA{i}{1} = 1; cphfA{i}{L(i)+1} = 1;
%     cphfy{i} = cell(L(i)+1,1);
%     cphfy{i}{1} = 1; cphfy{i}{L(i)+1} = 1;
% end;


last_sweep = false;

for swp=1:nswp
    % init check vector
%     rcchk = [1; checkrank*ones(d-1,1); 1];
%     rfchk = zeros(max(L)+1, d);
%     for i=1:d
%         rfchk(1:L(i)+1, i) = [1; checkrank*ones(L(i),1)];
%     end;

    dx_max = 0;
    res_max = 0;
    r_max = 0;
%     chk_res_max = 0;
    % bottom-to-top QR and phis
    for i=d:-1:1 % physical dims/core
        for j=1:L(i) % quantics dims
            cr = xf{i}{j};
            cr = reshape(cr, rfx(j,i)*n(j,i), rfx(j+1,i));
            [cr, rv] = qr(cr, 0);
            % What is our next core?
            if (j<L(i))
                % we are still on a "tooth"
                cr2 = xf{i}{j+1};
                cr2 = reshape(cr2, rfx(j+1,i), n(j+1,i)*rfx(j+2,i));
                cr2 = rv*cr2;
                rfx(j+1,i) = size(cr, 2);
                xf{i}{j} = reshape(cr, rfx(j,i), n(j,i), rfx(j+1,i));
                xf{i}{j+1} = reshape(cr2, rfx(j+1,i), n(j+1,i), rfx(j+2,i));
            else
                % We have to convlove rv to the tucker core
                cr2 = xc{i};
                cr2 = permute(cr2, [2, 1, 3]);
                cr2 = reshape(cr2, rfx(j+1,i), rcx(i)*rcx(i+1));
                cr2 = rv*cr2;
                rfx(j+1,i) = size(cr, 2);
                cr2 = reshape(cr2, rfx(j+1,i), rcx(i), rcx(i+1));
                cr2 = permute(cr2, [2, 1, 3]);
                xf{i}{j} = reshape(cr, rfx(j,i), n(j,i), rfx(j+1,i));
                xc{i} = cr2;
            end;
            % Update bottom phis
            cr = reshape(cr, rfx(j,i), n(j,i), rfx(j+1,i));
            if (j<L(i))
                phfA{i}{j+1} = compute_next_Phi(phfA{i}{j}, cr, Af{i}{j}, cr, 'lr');
                phfy{i}{j+1} = compute_next_Phi(phfy{i}{j}, cr, [], yf{i}{j}, 'lr');
            else
                phAfc{i} = compute_next_Phi(phfA{i}{j}, cr, Af{i}{j}, cr, 'lr');
                phyfc{i} = compute_next_Phi(phfy{i}{j}, cr, [], yf{i}{j}, 'lr');
            end;

            % check vector
%             ccr = ones(1, n(j,i), 1);
%             if (j<L(i))
%                 cphfA{i}{j+1} = compute_next_Phi(cphfA{i}{j}, ccr, Af{i}{j}, cr, 'lr');
%                 cphfy{i}{j+1} = compute_next_Phi(cphfy{i}{j}, ccr, [], yf{i}{j}, 'lr');
%             else
%                 cphAfc{i} = compute_next_Phi(cphfA{i}{j}, ccr, Af{i}{j}, cr, 'lr');
%                 cphyfc{i} = compute_next_Phi(cphfy{i}{j}, ccr, [], yf{i}{j}, 'lr');
%             end;
        end;
    end;

    % QRs and phis over the core
    % Project the system on the factors
    Acr = tt_matrix(Ac, Ac.n, ones(d,1));
%     cAcr = tt_matrix(Ac, Ac.n, ones(d,1));
    ycr = yc;
%     cycr = yc;
    for i=1:d
        Acr{i} = core_matrix(Ac{i}, phAfc{i});
        ycr{i} = core_vector(yc{i}, phyfc{i});
%         cAcr{i} = core_matrix(Ac{i}, cphAfc{i});
%         cycr{i} = core_vector(yc{i}, cphyfc{i});
    end;
    for i=d:-1:2
        rtx = rfx(L(i)+1, i);
        cr = xc{i};
        cr = reshape(cr, rcx(i), rtx*rcx(i+1));
        [cr, rv] = qr(cr.', 0);
        cr2 = xc{i-1};
        rtx2 = rfx(L(i-1)+1, i-1);
        cr2 =  reshape(cr2, rcx(i-1)*rtx2, rcx(i));
        cr2 = cr2*(rv.');
        rcx(i) = size(cr, 2);
        cr = reshape(cr.', rcx(i), rtx, rcx(i+1));
        xc{i-1} = reshape(cr2, rcx(i-1), rtx2, rcx(i));
        xc{i} = cr;
        % Update right phi
        phcA{i} = compute_next_Phi(phcA{i+1}, cr, Acr{i}, cr, 'rl');
        phcy{i} = compute_next_Phi(phcy{i+1}, cr, [], ycr{i}, 'rl');

        % check vector
%         ccr = 1;
%         % Update right phi
%         cphcA{i} = compute_next_Phi(cphcA{i+1}, ccr, cAcr{i}, cr, 'rl');
%         cphcy{i} = compute_next_Phi(cphcy{i+1}, ccr, [], cycr{i}, 'rl');
    end;

%     % DMRG over the core
%     % Compute the reduced matrix and rhs;
%     rta = zeros(d,1);
%     for i=1:d
%         rta(i) = rfA(L(i)+1,i);
%     end;
%     Acr = tt_matrix(Ac, rta, ones(d,1));
%     ycr = yc;
%     for i=1:d
%         a1 = Ac{i};
%         a1 = permute(a1, [2, 1, 3]);
%         a1 = reshape(a1, rfA(L(i)+1,i), rcA(i)*rcA(i+1));
%         ph2 = phfAb{i};
%         ph2 = permute(ph2, [1,3,2]);
%         ph2 = reshape(ph2, rfx(L(i)+1,i)*rfx(L(i)+1,i), rfA(L(i)+1,i));
%         a1 = ph2*a1;
%         a1 = reshape(a1, rfx(L(i)+1,i), rfx(L(i)+1,i), rcA(i), rcA(i+1));
%         a1 = permute(a1, [3, 1, 2, 4]);
%         Acr{i} = a1;
%
%         y1 = yc{i};
%         y1 = permute(y1, [2, 1, 3]);
%         y1 = reshape(y1, rfy(L(i)+1,i), rcy(i)*rcy(i+1));
%         ph2 = phfyb{i};
%         y1 = ph2*y1;
%         y1 = reshape(y1, rfx(L(i)+1,i), rcy(i), rcy(i+1));
%         y1 = permute(y1, [2, 1, 3]);
%         ycr{i} = y1;
%     end;
%     fprintf('=rake_solve========= core processing, sweep %d =======\n', swp);
%     xcold = xc;
%     xc = dmrg_solve2(Acr, ycr, tol, 'x0', xcold, 'max_full_size', max_full_size, 'nrestart', nrestart, 'gmres_iters', gmres_iters, 'nswp', 1);
%     res = norm(Acr*xcold-ycr)/norm(ycr);
%     dx = norm(xcold-xc)/norm(xc);
%     dx_max = max(dx_max, dx);
%     r_max = max([r_max; rank(xc)]);
%     fprintf('=rake_solve========= core res_prev: %3.3e, dx: %3.3e, rmax: %d\n', res, dx, max(rank(xc)));
%     res_max = max(res_max, res);

%     % Horisontal phis
%     rcx = xc.r;
%     for i=d:-1:2
%         % n here is in fact a tucker rank
%         cr2 = xc{i};
%         cr2 = reshape(cr2, rcx(i), rfx(L(i)+1,i)*rcx(i+1));
%         [cr2, rv]=qr(cr2.', 0);
%         cr3 = xc{i-1};
%         cr3 = reshape(cr3, rcx(i-1)*rfx(L(i-1)+1,i-1), rcx(i));
%         cr3 = cr3*(rv.');
%         rcx(i) = size(cr2, 2);
%         cr2 = cr2.'; % sizes rcx1, rtuck*rcx2
%         xc{i} = reshape(cr2, rcx(i), rfx(L(i)+1,i), rcx(i+1));
%         xc{i-1} = reshape(cr3, rcx(i-1), rfx(L(i-1)+1,i-1), rcx(i));
%
%         % Update right phis
%         cr2 = reshape(cr2, rcx(i), rfx(L(i)+1,i), rcx(i+1));
%         phcAr{i} = compute_next_Phi(phcAr{i+1}, cr2, Acr{i}, cr2, 'rl');
%         % New size: rcx1, rcA1, rcx1
%         phcyr{i} = compute_next_Phi(phcyr{i+1}, cr2, [], ycr{i}, 'rl');
%         % New size: rcx1, rcy1
%     end;

    % DMRG over the factors
    for i=1:d
        % Convolve the tucker block to the last physical
        % Preparing the matrix and rhs in this case smells like
        % manure.
        % phcAl - Ac{i} - phcAr     <- this has to be convolved
        %          |                <- this has to be convolved
        %         Af{i}{j}
        %          |
        %         Af{i}{j-1}
        %          |                <- this has to be convolved
        %         phfAb{j-1}
        rta = rfA(L(i)+1,i); rtx = rfx(L(i)+1,i); rty = rfy(L(i)+1,i); % rtchk = rfchk(L(i)+1,i);
        a1left = phcA{i};
        a1left = permute(a1left, [1,3,2]);
        a1left = reshape(a1left, rcx(i)*rcx(i), rcA(i));
        a1left = a1left*reshape(Ac{i}, rcA(i), rta*rcA(i+1));

%         % check
%         ca1left = cphcA{i};
%         ca1left = permute(ca1left, [1,3,2]);
%         ca1left = reshape(ca1left, rcchk(i)*rcx(i), rcA(i));
%         ca1left = ca1left*reshape(Ac{i}, rcA(i), rta*rcA(i+1));
        
        % old
%         % a1left is to use later to compute phcAl
%         ph2 = phcAr{i+1};
%         ph2 = permute(ph2, [2, 1, 3]);
%         ph2 = reshape(ph2, rcA(i+1), rcx(i+1)*rcx(i+1));
%         a1top = reshape(a1left, rcx(i)*rcx(i)*rta, rcA(i+1))*ph2;
%         a1top = reshape(a1top, rcx(i), rcx(i), rta, rcx(i+1), rcx(i+1));
%         a1top = permute(a1top, [1, 4, 2, 5, 3]);
%         % we need a1 as well, to check the  residual in tucker block
%         % splitting
%         a1top = reshape(a1top, rcx(i)*rcx(i+1)*rcx(i)*rcx(i+1), rta);
        % new
        a1left = reshape(a1left, rcx(i)*rcx(i), rta, rcA(i+1));
        a1left = permute(a1left, [2, 1, 3]);
        a1left = reshape(a1left, rta, rcx(i)*rcx(i)*rcA(i+1));

%         ca1left = reshape(ca1left, rcchk(i)*rcx(i), rta, rcA(i+1));
%         ca1left = permute(ca1left, [2, 1, 3]);
%         ca1left = reshape(ca1left, rta, rcchk(i)*rcx(i)*rcA(i+1));


        a2 = Af{i}{L(i)};
        a2 = reshape(a2, rfA(L(i),i)*n(L(i),i)*n(L(i),i), rta);
        % old
%         a2 = a2*a1top.';
%         a2 = reshape(a2, rfA(L(i),i), n(L(i),i), n(L(i),i), rcx(i)*rcx(i+1), rcx(i)*rcx(i+1));
%         a2 = permute(a2, [1, 2, 4, 3, 5]);
%         a2 = reshape(a2, rfA(L(i),i), n(L(i),i)*rcx(i)*rcx(i+1), n(L(i),i)*rcx(i)*rcx(i+1), 1);
        % new
        a2 = a2*a1left;
        a2 = reshape(a2, rfA(L(i),i), n(L(i),i), n(L(i),i), rcx(i), rcx(i), rcA(i+1));
        a2 = permute(a2, [1, 2, 4, 3, 5, 6]);
        a2 = reshape(a2, rfA(L(i),i), n(L(i),i)*rcx(i), n(L(i),i)*rcx(i), rcA(i+1));

%         % check
%         ca2 = Af{i}{L(i)};
%         ca2 = reshape(ca2, rfA(L(i),i)*n(L(i),i)*n(L(i),i), rta);
%         ca2 = ca2*ca1left;
%         ca2 = reshape(ca2, rfA(L(i),i), n(L(i),i), n(L(i),i), rcchk(i), rcx(i), rcA(i+1));
%         ca2 = permute(ca2, [1, 2, 4, 3, 5, 6]);
%         ca2 = reshape(ca2, rfA(L(i),i), n(L(i),i)*rcchk(i), n(L(i),i)*rcx(i), rcA(i+1));

        Afr = Af{i};
        Afr{L(i)} = a2;

%         cAfr = Af{i};
%         cAfr{L(i)} = ca2;

        y1left = phcy{i};
        y1left = y1left*reshape(yc{i}, rcy(i), rty*rcy(i+1));
        % y1left is to use later to compute phcyl
        ph2 = phcy{i+1};
        ph2 = ph2.';
        y1top = reshape(y1left, rcx(i)*rty, rcy(i+1))*ph2;
        y1top = reshape(y1top, rcx(i), rty, rcx(i+1));
        y1top = permute(y1top, [1, 3, 2]);
        y1top = reshape(y1top, rcx(i)*rcx(i+1), rty);
        y2 = yf{i}{L(i)};
        y2 = reshape(y2, rfy(L(i),i)*n(L(i),i), rty);
        y2 = y2*y1top.';
        y2 = reshape(y2, rfy(L(i),i), n(L(i),i)*rcx(i)*rcx(i+1), 1);
        yfr = yf{i};
        yfr{L(i)} = y2;

%         cy1left = cphcy{i};
%         cy1left = cy1left*reshape(yc{i}, rcy(i), rty*rcy(i+1));
%         % y1left is to use later to compute phcyl
%         ph2 = cphcy{i+1};
%         ph2 = ph2.';
%         cy1top = reshape(cy1left, rcchk(i)*rty, rcy(i+1))*ph2;
%         cy1top = reshape(cy1top, rcchk(i), rty, rcchk(i+1));
%         cy1top = permute(cy1top, [1, 3, 2]);
%         cy1top = reshape(cy1top, rcchk(i)*rcchk(i+1), rty);
%         cy2 = yf{i}{L(i)};
%         cy2 = reshape(cy2, rfy(L(i),i)*n(L(i),i), rty);
%         cy2 = cy2*cy1top.';
%         cy2 = reshape(cy2, rfy(L(i),i), n(L(i),i)*rcchk(i)*rcchk(i+1), 1);
%         cyfr = yf{i};
%         cyfr{L(i)} = cy2;

        x1 = xc{i};
        x1 = permute(x1, [2, 1, 3]);
        x1 = reshape(x1, rtx, rcx(i)*rcx(i+1));
        x2 = xf{i}{L(i)};
        x2 = reshape(x2, rfx(L(i),i)*n(L(i),i), rtx);
        x2 = x2*x1;
        x2 = reshape(x2, rfx(L(i),i), n(L(i),i)*rcx(i)*rcx(i+1), 1);
        xfr = xf{i};
        xfr{L(i)} = x2;

%          chk1 = chkc{i};
%          chk1 = permute(chk1, [2, 1, 3]);
%          chk1 = reshape(chk1, rtchk, rcchk(i)*rcchk(i+1));
%          chk2 = chkf{i}{L(i)};
%          chk2 = reshape(chk2, rfchk(L(i),i)*n(L(i),i), rtchk);
%          chk2 = chk2*chk1;
%          chk2 = reshape(chk2, rfchk(L(i),i), n(L(i),i)*rcchk(i)*rcchk(i+1), 1);
%          chkfr = chkf{i};
%          chkfr{L(i)} = chk2;

        curn = xfr.n;
%         curm = xfr.n;
%         curm(L(i)) = n(L(i),i); % *rcchk(i)*rcchk(i+1);
        curra = Afr.r;
        curry = yfr.r;
        currx = xfr.r;
%         currchk = [rfchk(1:L(i),i); 1];
%          currchk = ones(L(i)+1,1);
%          currchk = chkfr.r;
        % DMRG over the factor from L(i) to 1
        for j=L(i):-1:2
            % new
            if (j<L(i))
                Phi2 = phfA{i}{j+1};
%                 cPhi2 = cphfA{i}{j+1};
            else
                Phi2 = phcA{i+1};
%                 cPhi2 = cphcA{i+1};
            end;
            a2 = Afr{j};
            a1 = Afr{j-1};
            Phi1 = phfA{i}{j-1};
%             % check
%             ca2 = cAfr{j};
%             ca1 = cAfr{j-1};
%             cPhi1 = cphfA{i}{j-1};
            % old
%             a2 = phfAt{i}{j+1};
%             a2 = permute(a2, [2, 1, 3]);
%             a2 = reshape(a2, curra(j+1), currx(j+1)*currx(j+1));
%             a2 = reshape(Afr{j}, curra(j)*curn(j)*curn(j), curra(j+1))*a2;
%             a2 = reshape(a2, curra(j), curn(j), curn(j), currx(j+1), currx(j+1));
%             a2 = permute(a2, [2, 4, 3, 5, 1]);
%             a2 = reshape(a2, curn(j)*currx(j+1), curn(j)*currx(j+1), curra(j));
%             a1 = phfAb{i}{j-1};
%             a1 = permute(a1, [1, 3, 2]);
%             a1 = reshape(a1, currx(j-1)*currx(j-1), curra(j-1));
%             a1 = a1*reshape(Afr{j-1}, curra(j-1), curn(j-1)*curn(j-1)*curra(j));
%             a1 = reshape(a1, currx(j-1), currx(j-1), curn(j-1), curn(j-1), curra(j));
%             a1 = permute(a1, [1, 3, 2, 4, 5]);
%             a1 = reshape(a1, currx(j-1)*curn(j-1), currx(j-1)*curn(j-1), curra(j));

            y2 = phfy{i}{j+1}.';
            y2 = reshape(yfr{j}, curry(j)*curn(j), curry(j+1))*y2;
            y2 = reshape(y2, curry(j), curn(j)*currx(j+1));
            y2 = y2.';

            y1 = phfy{i}{j-1};
            y1 = y1*reshape(yfr{j-1}, curry(j-1), curn(j-1)*curry(j));
            y1 = reshape(y1, currx(j-1)*curn(j-1), curry(j));

            x2 = reshape(xfr{j}, currx(j), curn(j)*currx(j+1));
            x2 = x2.';
            x1 = reshape(xfr{j-1}, currx(j-1)*curn(j-1), currx(j));

%             % check
%             cy2 = cphfy{i}{j+1}.';
%             cy2 = reshape(cyfr{j}, curry(j)*curm(j), curry(j+1))*cy2;
%             cy2 = reshape(cy2, curry(j), curm(j)*currchk(j+1));
%             cy1 = cphfy{i}{j-1};
%             cy1 = cy1*reshape(cyfr{j-1}, curry(j-1), curm(j-1)*curry(j));
%             cy1 = reshape(cy1, currchk(j-1)*curm(j-1), curry(j));
%             cy = cy1*cy2;

            % new
            if (j==L(i))
                currx(j+1)=rcx(i+1);
%                 currchk(j+1)=rcchk(i+1);
                curn(j) = n(L(i),i)*rcx(i);
%                 curm(j) = n(L(i),i)*rcchk(i);
            end;

            if (verb>1)
                fprintf('=rake_solve2= swp %d, factor {%d}{%d}, ', swp, i, j);
            end;
            local_format = 'full';
            if (currx(j-1)*curn(j-1)*curn(j)*currx(j+1)>max_full_size2)
                local_format = 'tt';
            end;
            % old
%             [u,s,v,r,dx_max,res_max]=local_solve(a1, a2, y1, y2, x1, x2, ...
%                 currx(j-1), curn(j-1), curn(j), currx(j+1), curra(j), ...
%                 tol/sqrt(L(i))/sqrt(d)/2, res_max, dx_max, ...
%                 local_format, max_full_size, nrestart, gmres_iters, verb);
            % new
            [u,s,v,r,dx_max,res_max]=local_solve(Phi1,a1, a2, Phi2, y1, y2, x1, x2, ...
                currx(j-1), curn(j-1), curn(j), currx(j+1), curra(j), ...
                tol2/sqrt(sum(L)), res_max, dx_max, resid_damp_loc, trunc_norm, ...
                local_format, max_full_size, nrestart, gmres_iters, verb);
            % old rounding tol: tol2/sqrt(L(i))/sqrt(d)/2
	    r = min(rmax, r);
	    u = u(:,1:r); s = s(1:r,1:r); v = v(:,1:r);
            u = u*s;

%             % check
%             Asol = bfun3(cPhi1,ca1,ca2,cPhi2, u*(v.'));
%             chk_res_max = max(chk_res_max, norm(Asol-cy(:))/norm(cy(:)));

            % kick
            if (~last_sweep)
%                 Axprev = bfun3(Phi1,a1, a2, Phi2, x1*x2.');
%                 Axprev = reshape(Axprev, currx(j-1)*curn(j-1), curn(j)*currx(j+1));
%                 [unew,snew,vnew]=svd(Axprev, 'econ');
%                 v = reort(v, vnew(:,1:min(kickrank, size(vnew,2))));
%                 v = reort(v, rand(curn(j)*currx(j+1), kickrank));
                [v,rv]=qr([v, rand(curn(j)*currx(j+1), kickrank)], 0);
                radd = kickrank;
                u = [u, zeros(currx(j-1)*curn(j-1), radd)];
                u = u*(rv.');
            end;
            r = size(v,2);
            xfr{j} = reshape(v.', r, curn(j), currx(j+1));
            xfr{j-1} = reshape(u, currx(j-1), curn(j-1), r);
            currx(j) = r;
            rfx(j,i)=r;
            r_max = max(r_max, r);
            % old
            % new phis
%             a2 = permute(a2, [1, 3, 2]);
%             a2 = reshape(a2, curn(j)*currx(j+1)*curra(j), curn(j)*currx(j+1));
%             a2 = a2*v;
%             a2 = reshape(a2, curn(j)*currx(j+1), curra(j)*r);
%             a2 = (v')*a2;
%             phfAt{i}{j} = reshape(a2, r, curra(j), r);
            phfy{i}{j} = (v')*y2;
            %new
            phfA{i}{j} = compute_next_Phi(Phi2, xfr{j}, a2, xfr{j}, 'rl');

%             % check vector
%             ccr = ones(n(j,i), 1);
%             cphfy{i}{j} = (ccr')*(cy2.');
%             cphfA{i}{j} = compute_next_Phi(cPhi2, ccr.', ca2, xfr{j}, 'rl');
        end;

%         fprintf('=rake_solve========= factor {%d} processing, Sweep %d =======\n', i, swp);
%         xfrold = xfr;
%         xfr = dmrg_solve2(Afr, yfr, tol, 'x0', xfrold, 'max_full_size', max_full_size, 'nrestart', nrestart, 'gmres_iters', gmres_iters, 'nswp', 1);
%         res = norm(Afr*xfrold-yfr)/norm(yfr);
%         dx = norm(xfrold-xfr)/norm(xfr);
%         dx_max = max(dx_max, dx);
%         r_max = max([r_max; rank(xfr)]);
%         fprintf('=rake_solve========= factor {%d} res_prev: %3.3e, dx: %3.3e, rmax: %d\n', i, res, dx, max(rank(xfr)));
%         res_max = max(res_max, res);
%         rfx(1:L(i)+1,i) = xfr.r;

        % We have to split the tucker block, and compute new phf*b
        for j=1:L(i)
            cr = xfr{j};
            % What is our next core?
            if (j<L(i))
                % we are still on a "tooth"
                cr = reshape(cr, rfx(j,i)*n(j,i), rfx(j+1,i));
                [cr, rv] = qr(cr, 0);
                cr2 = xfr{j+1};
                ncur = size(cr2, 2);
                r3cur = size(cr2, 3);
                cr2 = reshape(cr2, rfx(j+1,i), ncur*r3cur);
                cr2 = rv*cr2;
                rfx(j+1,i) = size(cr, 2);
                xfr{j} = reshape(cr, rfx(j,i), n(j,i), rfx(j+1,i));
                xfr{j+1} = reshape(cr2, rfx(j+1,i), ncur, r3cur);
            else
                % We have to split the tucker core and update the tucker
                % rank
                cr = reshape(cr, rfx(j,i)*n(j,i), rcx(i)*rcx(i+1));
                [u,s,v]=svd(cr, 'econ');
                if (strcmp(trunc_norm, 'fro'))
                    r = my_chop2(diag(s), tol/sqrt(L(i))/sqrt(d)/2*norm(diag(s)));
                else
                % Prepare the local matrix and rhs for residue check
                % new
                Phi2 = phcA{i+1};
                curA2 = reshape(a1left, rta, rcx(i), rcx(i), rcA(i+1));
                curA1 = Af{i}{L(i)};
                Phi1 = phfA{i}{j};
                % old
%                 curA = cell(2,1);
%                 curA{2} = reshape(a1top, rcx(i)*rcx(i+1), rcx(i)*rcx(i+1), rta);
%                 curA{1} = Af{i}{L(i)};
%                 curA{1} = reshape(curA{1}, rfA(j,i), n(j,i)*n(j,i)*rta);
%                 ph2 = phfAb{i}{j};
%                 ph2 = permute(ph2, [1,3,2]);
%                 ph2 = reshape(ph2, rfx(j,i)*rfx(j,i), rfA(j,i));
%                 curA{1} = ph2*curA{1};
%                 curA{1} = reshape(curA{1}, rfx(j,i), rfx(j,i), n(j,i), n(j,i), rta);
%                 curA{1} = permute(curA{1}, [1, 3, 2, 4, 5]);
%                 curA{1} = reshape(curA{1}, rfx(j,i)*n(j,i), rfx(j,i)*n(j,i), rta);

                rhs = reshape(yf{i}{j}, rfy(j,i), n(j,i)*rty);
                rhs = phfy{i}{j}*rhs;
                rhs = reshape(rhs, rfx(j,i)*n(j,i), rty);
                rhs = rhs*(y1top.');
                rhs = reshape(rhs, rfx(j,i)*n(j,i)*rcx(i)*rcx(i+1),1);
                r = 1;
                normy = norm(rhs);
                % old
%                 res_true = norm(bfun2(curA, cr, rfx(j,i), n(j,i), rcx(i), rcx(i+1), rfx(j,i), n(j,i), rcx(i), rcx(i+1))-rhs)/normy;
                % new
                res_true = norm(bfun3(Phi1, curA1, curA2, Phi2, cr)-rhs)/normy;
                while (r<=size(s,1))
                    cursol = u(:,1:r)*s(1:r,1:r)*(v(:,1:r)');
                    % old
%                     res = norm(bfun2(curA, cursol, rfx(j,i), n(j,i), rcx(i), rcx(i+1), rfx(j,i), n(j,i), rcx(i), rcx(i+1))-rhs)/normy;
                    % new
                    res = norm(bfun3(Phi1, curA1, curA2, Phi2, cursol)-rhs)/normy;
                    if (res<max(tol/sqrt(sum(L)), res_true*resid_damp_loc))
                        % old: tol/sqrt(L(i))/sqrt(d)/2
                        break;
                    end;
                    r = r+1;
                end;
                end;
                
                if (verb>1)
                    fprintf('=rake_solve2= swp %d, tuckerrank {%d}, res: %3.3e, r: %d\n', swp, i, res, r);
                end;
                r = min(rmax, r);
                u = u(:,1:r);
                v = conj(v(:,1:r));
                s = s(1:r,1:r);
                v = v*s;
                if (~last_sweep)
%                     Axprev = bfun3(Phi1, curA1, curA2, Phi2, cr);
%                     Axprev = reshape(Axprev, rfx(j,i)*n(j,i), rcx(i)*rcx(i+1));
%                     [unew,snew,vnew]=svd(Axprev, 'econ');
%                     u = reort(u, unew(:,1:min(kickrank, size(unew,2))));                
%                     u = reort(u, rand(rfx(j,i)*n(j,i), kickrank));
                    [u,rv]=qr([u, rand(rfx(j,i)*n(j,i), kickrank)], 0);
                    radd = kickrank;
                    v = [v, zeros(rcx(i)*rcx(i+1), radd)];
                    v = v*(rv.');
                end;
                r = size(u,2);
                rfx(j+1,i) = r;
                cr = u;
                xfr{j} = reshape(cr, rfx(j,i), n(j,i), r);
                v = reshape(v.', r, rcx(i), rcx(i+1));
                v = permute(v, [2, 1, 3]);
                xc{i} = v;
                r_max = max(r_max, r);
            end;
            % Update bottom phis
            cr = reshape(cr, rfx(j,i), n(j,i), rfx(j+1,i));
            if (j<L(i))
                phfA{i}{j+1} = compute_next_Phi(phfA{i}{j}, cr, Af{i}{j}, cr, 'lr');
                phfy{i}{j+1} = compute_next_Phi(phfy{i}{j}, cr, [], yf{i}{j}, 'lr');
            else
                phAfc{i} = compute_next_Phi(phfA{i}{j}, cr, Af{i}{j}, cr, 'lr');
                phyfc{i} = compute_next_Phi(phfy{i}{j}, cr, [], yf{i}{j}, 'lr');
            end;

%             % check vector
%             ccr = ones(1, n(j,i), 1);
%             if (j<L(i))
%                 cphfA{i}{j+1} = compute_next_Phi(cphfA{i}{j}, ccr, Af{i}{j}, cr, 'lr');
%                 cphfy{i}{j+1} = compute_next_Phi(cphfy{i}{j}, ccr, [], yf{i}{j}, 'lr');
%             else
%                 cphAfc{i} = compute_next_Phi(cphfA{i}{j}, ccr, Af{i}{j}, cr, 'lr');
%                 cphyfc{i} = compute_next_Phi(cphfy{i}{j}, ccr, [], yf{i}{j}, 'lr');
%             end;
        end;
        xf{i} = xfr;

        % Now, perform the DMRG over the core to the next factor, update phc*l as well
        if (i<d)
            rtx = rfx(L(i)+1,i); rtx2 = rfx(L(i+1)+1,i+1);
            rx1 = rcx(i); rx2 = rcx(i+1); rx3 = rcx(i+2);
            ra2 = rcA(i+1); ra3 = rcA(i+2);
            ry2 = rcy(i+1); ry3 = rcy(i+2);

%             rtchk = rfchk(L(i)+1,i); rtchk2 = rfchk(L(i+1)+1,i+1);
%             rchk1 = rcchk(i); rchk3 = rcchk(i+2);

            % Acr and ycr are okey, except the i-th core
            % new
            Phi1 = phcA{i};
            a1 = core_matrix(Ac{i}, phAfc{i});
            a2 = Acr{i+1};
            Phi2 = phcA{i+2};

%             cPhi1 = cphcA{i};
%             ca1 = core_matrix(Ac{i}, cphAfc{i});
%             ca2 = cAcr{i+1};
%             cPhi2 = cphcA{i+2};
            % old
%             a1 = reshape(a1left, rx1*rx1, rta, ra2);
%             a1 = core_matrix(a1, phAfc{i});
%             a1 = reshape(a1, rx1, rx1, rtx, rtx, ra2);
%             a1 = permute(a1, [1, 3, 2, 4, 5]);
%             a1 = reshape(a1, rx1*rtx, rx1*rtx, ra2);
%             a2 = phcAr{i+2};
%             a2 = permute(a2, [2, 1, 3]);
%             a2 = reshape(a2, ra3, rx3*rx3);
%             a2 = reshape(Acr{i+1}, ra2*rtx2*rtx2, ra3)*a2;
%             a2 = reshape(a2, ra2, rtx2, rtx2, rx3, rx3);
%             a2 = permute(a2, [2, 4, 3, 5, 1]);
%             a2 = reshape(a2, rtx2*rx3, rtx2*rx3, ra2);

            y1 = reshape(y1left, rx1, rty, ry2);
            y1 = core_vector(y1, phyfc{i});
            y1 = reshape(y1, rx1*rtx, ry2);

            y2 = phcy{i+2}.';
            y2 = reshape(ycr{i+1}, ry2*rtx2, ry3)*y2;
            y2 = reshape(y2, ry2, rtx2*rx3);
            y2 = y2.';

            x1 = reshape(xc{i}, rx1*rtx, rx2);
            x2 = reshape(xc{i+1}, rx2, rtx2*rx3);
            x2 = x2.';

            if (verb>1)
                fprintf('=rake_solve2= swp %d, core {%d}, ', swp, i);
            end;
            local_format = 'full';
            if (rx1*rtx*rtx2*rx3>max_full_size2)
                local_format = 'tt';
            end;
            % new
            [u,s,v,r,dx_max,res_max]=local_solve(Phi1, a1, a2, Phi2, y1, y2, x1, x2, ...
                rx1, rtx, rtx2, rx3, ra2, ...
                tol2/sqrt(sum(L)), res_max, dx_max, resid_damp_loc, trunc_norm, ...
                local_format, max_full_size, nrestart, gmres_iters, verb);
            % old tol: tol2/sqrt(d)/2
            % old
%             [u,s,v,r,dx_max,res_max]=local_solve(a1, a2, y1, y2, x1, x2, ...
%                 rx1, rtx, rtx2, rx3, ra2, ...
%                 tol/sqrt(d)/2, res_max, dx_max, ...
%                 local_format, max_full_size, nrestart, gmres_iters, verb);
            r = min(rmax, r);
            u = u(:,1:r); s = s(1:r,1:r); v = v(:,1:r);
            v = v*s;

%             % check
%             Asol = bfun3(cPhi1, ca1, ca2, cPhi2, u*(v.'));
%             cy1 = reshape(cy1left, rchk1, rty, ry2);
%             cy1 = core_vector(cy1, cphyfc{i});
%             cy1 = reshape(cy1, rchk1*rtchk, ry2);
%             cy2 = cphcy{i+2}.';
%             cy2 = reshape(cycr{i+1}, ry2*rtchk2, ry3)*cy2;
%             cy2 = reshape(cy2, ry2, rtchk2*rchk3);
%             cy = cy1*cy2;
%             chk_res_max = max(chk_res_max, norm(Asol-cy(:))/norm(cy(:)));

            % kick
            if (~last_sweep)
%                 Axprev = bfun3(Phi1, a1, a2, Phi2, x1*x2.');
%                 Axprev = reshape(Axprev, rx1*rtx, rtx2*rx3);
%                 [unew,snew,vnew]=svd(Axprev, 'econ');    
%                 u = reort(u, unew(:,1:min(kickrank, size(unew,2))));                
%                 u = reort(u, rand(rx1*rtx, kickrank));
                [u,rv]=qr([u, rand(rx1*rtx, kickrank)], 0);
                radd = kickrank;
                v = [v, zeros(rtx2*rx3, radd)];
                v = v*(rv.');
            end;
            r = size(u,2);
            xc{i} = reshape(u, rx1, rtx, r);
            xc{i+1} = reshape(v.', r, rtx2, rx3);
            rcx(i+1)=r;
            r_max = max(r_max, r);
            % new phis
            % old
%             a1 = permute(a1, [1, 3, 2]);
%             a1 = reshape(a1, rx1*rtx*ra2, rx1*rtx);
%             a1 = a1*u;
%             a1 = reshape(a1, rx1*rtx, ra2*r);
%             a1 = (u')*a1;
%             phcAl{i+1} = reshape(a1, r, ra2, r);
            phcy{i+1} = (u')*y1;
            % new
            phcA{i+1} = compute_next_Phi(phcA{i}, xc{i}, a1, xc{i}, 'lr');

%             ccr = 1;
%             cphcy{i+1} = (ccr')*cy1;
%             ccr = reshape(ccr, rchk1, rtchk, rcchk(i+1));
%             cphcA{i+1} = compute_next_Phi(cphcA{i}, ccr, ca1, xc{i}, 'lr');


%             cr = xc{i};
%             cr = reshape(cr, rcx(i)*rtx, rcx(i+1));
%             [cr, rv]=qr(cr, 0);
%             cr2 = xc{i+1};
%             cr2 = reshape(cr2, rcx(i+1), rfx(L(i+1)+1,i+1)*rcx(i+2));
%             cr2 = rv*cr2;
%             rcx(i+1) = size(cr, 2);
%             xc{i+1} = reshape(cr2, rcx(i+1), rfx(L(i+1)+1,i+1), rcx(i+2));
%             xc{i} = reshape(cr, rcx(i), rtx, rcx(i+1));
%
%             % We have a1left, y1left of sizes rcx(i)(^2), rtuck, rc*(i+1)
%             curph = reshape(a1left, rcx(i)*rcx(i), rta, rcA(i+1));
%             curph = permute(curph, [2, 1, 3]);
%             curph = reshape(curph, rta, rcx(i)*rcx(i)*rcA(i+1));
%             ph2 = phfAb{i};
%             ph2 = permute(ph2, [1, 3, 2]);
%             ph2 = reshape(ph2, rtx*rtx, rta);
%             curph = ph2*curph;
%             curph = reshape(curph, rtx, rtx, rcx(i), rcx(i), rcA(i+1));
%             curph = permute(curph, [3, 1, 5, 4, 2]);
%             curph = reshape(curph, rcx(i)*rtx*rcA(i+1), rcx(i)*rtx);
%             curph = curph*cr;
%             curph = reshape(curph, rcx(i)*rtx, rcA(i+1)*rcx(i+1));
%             curph = (cr')*curph;
%             phcAl{i+1} = reshape(curph, rcx(i+1), rcA(i+1), rcx(i+1));
%
%             curph = reshape(y1left, rcx(i), rty, rcy(i+1));
%             curph = permute(curph, [2, 1, 3]);
%             curph = reshape(curph, rty, rcx(i)*rcy(i+1));
%             ph2 = phfyb{i};
%             curph = ph2*curph;
%             curph = reshape(curph, rtx, rcx(i), rcy(i+1));
%             curph = permute(curph, [2, 1, 3]);
%             curph = reshape(curph, rcx(i)*rtx, rcy(i+1));
%             curph = (cr')*curph;
%             phcyl{i+1} = curph;
        end;
    end;

    if (verb>0)
        real_res = NaN;
%  	x = qtt_tucker;
%  	x.dphys = d;
%  	x.tuck = xf;
%  	x.core = xc;
%
%  	real_res = mvrk(A, x, tol, 'verb', 0);
%  %          real_res = A*x;
%  	real_res = norm(real_res-y)/norm(y);
        fprintf('\n=rake_solve2= swp %d, dx_max: %3.3e, res_max: %3.3e, r_max: %d, real_res: %3.3e\n\n', swp, dx_max, res_max, r_max, real_res);
    end;
    if (last_sweep)
        break;
    end;
    if (strcmp(trunc_norm, 'fro'))
        if (dx_max<tol)||(swp==nswp-1)
            last_sweep = true;
        end;
    else
        if (res_max<tol)||(swp==nswp-1)
            last_sweep=true;
        end;
%         if (chk_res_max<tol)||(swp==nswp-1)
%             last_sweep = true;
%         end;
    end;
end;

x = qtt_tucker;
x.dphys = d;
x.tuck = xf;
x.core = xc;
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

function [A] = core_matrix(core_block, Phi_factor)
% Computes the matrix to use in DMRG by the core, by convolving the core of
% the Tucker block and the last factor Phi matrix over the Tucker rank
% core_block is rc1 x rtuck x rc2, and
% Phi_factor is n1 x rtuck x n2
% A is then rc1 x n1 x n2 x rc2

r1 = size(core_block, 1); rtuck = size(core_block, 2); r2 = size(core_block, 3);
n1 = size(Phi_factor, 1); n2 = size(Phi_factor, 3);

Phi_factor = permute(Phi_factor, [1, 3, 2]);
Phi_factor = reshape(Phi_factor, n1*n2, rtuck);
A = permute(core_block, [2, 1, 3]);
A = reshape(A, rtuck, r1*r2);
A = Phi_factor*A;
A = reshape(A, n1, n2, r1, r2);
A = permute(A, [3, 1, 2, 4]);

end


function [A] = core_vector(core_block, Phi_factor)
% Computes the vector to use in DMRG by the core, by convolving the core of
% the Tucker block and the last factor Phi matrix over the Tucker rank
% core_block is rc1 x rtuck x rc2, and
% Phi_factor is n1 x rtuck
% A is then rc1 x n1 x rc2

r1 = size(core_block, 1); rtuck = size(core_block, 2); r2 = size(core_block, 3);
n1 = size(Phi_factor, 1);

A = permute(core_block, [2, 1, 3]);
A = reshape(A, rtuck, r1*r2);
A = Phi_factor*A;
A = reshape(A, n1, r1, r2);
A = permute(A, [2, 1, 3]);

end


function [y]=bfun2(B, x, rxm1, m1, m2, rxm3, rxn1, k1, k2, rxn3)
% Computes (B{1} \otimes B{2})x
% B{1} is of sizes rxn1*k1, rxm1*m1, rB
% B{2} is of sizes k2*rxn3, m2*rxm3, rB
rB=size(B{1},3);
x = reshape(x, rxm1*m1, m2*rxm3);
B1 = permute(B{1}, [3 1 2]);
B1 = reshape(B1, rB*rxn1*k1, rxm1*m1);
y = B1*x; % size rB*rxn1*k1,m2*rxm3  % cplx rB*rx^3*n^3
y = reshape(y, rB, rxn1*k1, m2*rxm3);
y = permute(y, [3 1 2]);
y = reshape(y, m2*rxm3*rB, rxn1*k1);
B2 = reshape(B{2}, k2*rxn3, m2*rxm3*rB);
y = B2*y; % size k2*rxn3,rxn1*k1 % cplx rB*rx^3*n^3
y = reshape(y.', rxn1*k1*k2*rxn3, 1);
end


function [y]=bfun3(Phi1,B1,B2,Phi2, x)
% Computes (Phi1 * B1 * B2 * Phi2)*x
% Phi1 is of sizes ry1, rB1, rx1
% B1 is of sizes rB1, k1, m1, rB2
% B2 is of sizes rB2, k2, m2, rB3
% Phi2 is of sizes ry3, rB3, rx3
ry1 = size(Phi1,1); ry3 = size(Phi2,1);
rx1 = size(Phi1,3); rx3 = size(Phi2,3);
rb1=size(B1,1); rb2=size(B1,4); rb3 = size(B2, 4);
m1 = size(B1,3); m2 = size(B2,3);
k1 = size(B1,2); k2 = size(B2,2);

y = reshape(x, rx1, m1*m2*rx3);
Phi1 = reshape(Phi1, ry1*rb1, rx1);
y = Phi1*y; % size ry1*rb1,m1*m2*rx3 % cplx rb*rx^3*m^2
y = reshape(y, ry1, rb1*m1, m2, rx3);
y = permute(y, [2, 1, 3, 4]);
y = reshape(y, rb1*m1, ry1*m2*rx3);
B1 = permute(B1, [2, 4, 1, 3]);
B1 = reshape(B1, k1*rb2, rb1*m1);
y = B1*y; % size k1*rb2, ry1*m2*rx3 % cplx rb^2*rx^2*n^3
y = reshape(y, k1, rb2, ry1, m2, rx3);
y = permute(y, [2, 4, 3, 1, 5]);
y = reshape(y, rb2*m2, ry1*k1*rx3);
B2 = permute(B2, [2, 4, 1, 3]);
B2 = reshape(B2, k2*rb3, rb2*m2);
y = B2*y; % size k2*rb3, ry1*k1*rx3 % cplx rb^2*rx^2*n^3
y = reshape(y, k2, rb3, ry1*k1, rx3);
y = permute(y, [2, 4, 3, 1]);
y = reshape(y, rb3*rx3, ry1*k1*k2);
Phi2 = reshape(Phi2, ry3, rb3*rx3);
y = Phi2*y; % size ry3, ry1*k1*k2 % cplx rb*rx^3*n^2
y = y.';
y = reshape(y, ry1*k1*k2*ry3, 1);
end


% old
% function [u,s,v,r,dx_max,res_max]=local_solve(a1, a2, y1, y2, x1, x2, ...
% new
function [u,s,v,r,dx_max,res_max]=local_solve(Phi1,a1, a2, Phi2, y1, y2, x1, x2, ...
    rx1, n1, n2, rx3, ra2, ...
    real_tol, res_max, dx_max,  resid_damp, trunc_norm, ...
    local_format, max_full_size, nrestart, gmres_iters, verb)

if (strcmp(local_format, 'full'))
    sol_prev = x1*(x2.');
    sol_prev = sol_prev(:);
    rhs = y1*(y2.');
    rhs = rhs(:);
    normy = norm(rhs);
    if (rx1*n1*n2*rx3<max_full_size)
        % new
        B = permute(Phi1, [1,3,2]);
        B = reshape(B, rx1*rx1, size(a1,1));
        a1 = reshape(a1, size(a1,1), n1*n1*ra2);
        B = B*a1;
        B = reshape(B, rx1, rx1, n1, n1, ra2);
        B = permute(B, [1, 3, 2, 4, 5]);
        B = reshape(B, rx1*n1*rx1*n1, ra2);
        ra3 = size(a2,4);
        a2 = reshape(a2, ra2, n2*n2*ra3);
        B = B*a2;
        B = reshape(B, rx1*n1*rx1*n1*n2*n2, ra3);
        Phi2 = permute(Phi2, [2, 1, 3]);
        Phi2 = reshape(Phi2, ra3, rx3*rx3);
        B = B*Phi2;
        B = reshape(B, rx1*n1, rx1*n1, n2, n2, rx3, rx3);
        B = permute(B, [1, 3, 5, 2, 4, 6]);
        % old
%         B = reshape(a1, rx1*n1*rx1*n1, ra2);
%         B = B*(reshape(a2, n2*rx3*n2*rx3, ra2).');
%         B = reshape(B, rx1*n1,  rx1*n1, n2*rx3, n2*rx3);
%         B = permute(B, [1, 3, 2, 4]);

        B = reshape(B, rx1*n1*n2*rx3, rx1*n1*n2*rx3);
%         nB = norm(B, 'fro');
%         B = B/nB;
%         rhs = rhs/nB;
%          sol = (B'*B+1e-15*eye(size(B))) \ (B'*rhs);
%          
%          B = B*nB;
%          rhs = rhs*nB;
         sol = B \ rhs;
               
        
%          sol = (B'*B + 1e-16*norm(B'*B, 'fro')) \ (B'*rhs);
%         sol = pinv(B)*rhs;
        flg = 0;
%         [sol,flg]=gmres(B, rhs, ...
%             nrestart, real_tol/resid_damp, gmres_iters, [], [], sol);

        res_true = norm(B*sol-rhs)/normy;
        res_prev = norm(B*sol_prev-rhs)/normy;
        
%         if (res_true/res_prev>1)&&(res_true>real_tol)
% %             If the direct solution sucked
%             [sol,flg]=gmres(B, rhs, ...
%                 nrestart, real_tol/resid_damp, gmres_iters, [], [], sol_prev);
% %             sol = sol_prev;
%             res_true = norm(B*sol-rhs)/normy;
%         else
%             flg = 0;
%         end;
        
        if (flg>0)
            fprintf('--warn-- gmres did not converge\n');
        end;         
    else
        B = cell(2,1);
        B{1} = a1;
        B{2} = a2;
        % old
%         res_prev = norm(bfun2(B, sol_prev, rx1, n1, n2, rx3, rx1, n1, n2, rx3)-rhs)/normy;
        % new
        drhs = bfun3(Phi1, a1, a2, Phi2, sol_prev)-rhs;
        res_prev = norm(drhs)/normy;

%         sol = als_solve_rx_2(B, rhs, real_tol, 10, sol_prev, [], 3);
%         [sol,flg]=bicgstab(@(v)bfun2(B, v, rx1, n1, n2, rx3, rx1, n1, n2, rx3), rhs, ...
%             max(real_tol,res_prev*0.1), nrestart*gmres_iters, [], [], sol_prev);
        % old
%         [sol,flg]=gmres(@(v)bfun2(B, v, rx1, n1, n2, rx3, rx1, n1, n2, rx3), rhs, ...
%             nrestart, max(real_tol,res_prev*0.05), gmres_iters, [], [], sol_prev);
        % new
        [dsol,flg]=gmres(@(v)bfun3(Phi1, a1, a2, Phi2, v), drhs, ...
            nrestart, min(real_tol/resid_damp/res_prev,1), gmres_iters);
        sol = sol_prev-dsol;
        if (flg>0)
            fprintf('--warn-- gmres did not converge\n');
        end;
        % old
%         res_true = norm(bfun2(B, sol, rx1, n1, n2, rx3, rx1, n1, n2, rx3)-rhs)/normy;
        % new
        res_true = norm(bfun3(Phi1, a1, a2, Phi2, sol)-rhs)/normy;
    end;

    if ((res_prev/res_true)<resid_damp)&&(res_true>real_tol/resid_damp)
        fprintf('--warn-- the residual damp by gmres was smaller than in the truncation\n');
    end;

    dx = norm(sol-sol_prev)/norm(sol);
    dx_max = max(dx_max, dx);
    if (rx1*n1*n2*rx3<max_full_size)
        Bx = B*sol; Bx_prev = B*sol_prev;
    else
        Bx = bfun3(Phi1, a1, a2, Phi2, sol);
        Bx_prev = bfun3(Phi1, a1, a2, Phi2, sol_prev);
    end;
    res_max = max(res_max, norm(Bx-Bx_prev)/norm(Bx));
%      res_max = max(res_max, res_prev);

    sol = reshape(sol, rx1*n1, n2*rx3);
    [u,s,v]=svd(sol, 'econ');
    s = diag(s);
    if (strcmp(trunc_norm, 'fro'))
        r = my_chop2(s, norm(s)*max(real_tol, res_true*resid_damp));
    else
        r1 = 1; r2 = numel(s); r = round((r1+r2)/2);
        while (r2-r1>1)
            cursol = u(:,1:r)*diag(s(1:r))*(v(:,1:r)');
            if (rx1*n1*n2*rx3<max_full_size)
                res = norm(B*cursol(:)-rhs)/normy;
            else
                % old
                %             res = norm(bfun2(B, cursol, rx1, n1, n2, rx3, rx1, n1, n2, rx3)-rhs)/normy;
                % new
                res = norm(bfun3(Phi1, a1, a2, Phi2, cursol)-rhs)/normy;
            end;
            if (res<max(real_tol, res_true*resid_damp))
                r2 = r;
            else
                r1 = r;
            end;
            r = round((r1+r2)/2);
        end;
        %     r = 1;
        while (r<=numel(s))
            cursol = u(:,1:r)*diag(s(1:r))*(v(:,1:r)');
            if (rx1*n1*n2*rx3<max_full_size)
                res = norm(B*cursol(:)-rhs)/normy;
            else
                % old
                %             res = norm(bfun2(B, cursol, rx1, n1, n2, rx3, rx1, n1, n2, rx3)-rhs)/normy;
                % new
                res = norm(bfun3(Phi1, a1, a2, Phi2, cursol)-rhs)/normy;
            end;
            if (res<max(real_tol, res_true*resid_damp))
                break;
            end;
            r = r+1;
        end;
    end;
    r = min(r, numel(s));
    if (verb>1)
        fprintf('dx: %3.3e, res: %3.3e, res_prev: %3.3e, r: %d\n', dx, res, res_prev, r);
    end;

    s = diag(s(1:r));
    u = u(:,1:r);
    v = conj(v(:,1:r));
else
    % implement tt-gmres here
    B = cell(2,1);
    B{1} = a1;
    B{2} = a2;
%     iB = tt_minres_selfprec(B, 1e-1, 1e-2, 10, 'right');
    iB = [];
    sol_prev = cell(2,1);
    sol_prev{1} = x1;
    sol_prev{2} = x2;
    rhs = cell(2,1);
    rhs{1} = y1;
    rhs{2} = y2;
    normy = tt_dist3(rhs, tt_scal(rhs,0));
    drhs = tt_mv(B, sol_prev);
    res_prev = tt_dist3(drhs, rhs)/normy;
    drhs = tt_add(rhs, tt_scal(drhs, -1));
    drhs = tt_compr2(drhs, real_tol);
    dsol = tt_gmres(B, drhs, real_tol*resid_damp_loc/res_prev, gmres_iters, nrestart, real_tol, real_tol, iB, [], [], [], 1);
%     sol = tt_gmres(B, rhs, real_tol*2, gmres_iters*5, nrestart/5, real_tol, real_tol, iB, [], [], sol_prev, 1);
    sol = tt_add(sol_prev, dsol);
    sol = tt_compr2(sol, real_tol);
    normsol = tt_dist3(sol, tt_scal(sol,0));
    dx = tt_dist3(sol, sol_prev)/normsol;
    res = tt_dist3(tt_mv(B, sol), rhs)/normy;

    dx_max = max(dx_max, dx);
    res_max = max(res_max, res_prev);
    [v, s]=qr(sol{2}, 0);
    sol{1} = sol{1}*(s.');
    [u, s]=qr(sol{1}, 0);
    r = size(sol{1},2);
    if (verb>1)
        fprintf('dx: %3.3e, res: %3.3e, res_prev: %3.3e, r: %d\n', dx, res, res_prev, r);
    end;
end;
end
