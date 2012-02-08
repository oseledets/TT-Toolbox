function [y]=mvrk(A, x, eps, varargin)
% Computes matvec in the QTT-Tucker "rake" format 
%   [Y]=MVRK(A, X, EPS, OPTIONS)
% Computes the MatVec y=Ax in the qtt_tucker "rake" format using the DMRG
% approach. Options are provided in form
%   'PropertyName1',PropertyValue1,'PropertyName2',PropertyValue2 and so
%   on. The parameters are set to default (in brackets in the following) 
%   The list of option names and default values is:
%       o kickrank -- the additional ranks, the larger the more robust the
%       method is, but the complexity increases [2]
%       o nswp - maximal number of DMRG sweeps [20]
%       o verb - verbosity level, 0-silent, 1-sweep info, 2-block info [1]
%       o y0 - initial approximation to the product [random-rank 2
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
kickrank = 2;
verb = 1;
y = [];

for i=1:2:length(varargin)-1
    switch lower(varargin{i})
        case 'nswp'
            nswp=varargin{i+1};
%         case 'rmax'
%             rmax=lower(varargin{i+1});
        case 'y0'
            y=varargin{i+1};
        case 'verb'
            verb=varargin{i+1};
        case 'kickrank'
            kickrank=varargin{i+1};
    end;
end;

d = x.dphys;
xc = x.core;
xf = x.tuck;
Af = A.tuck;
Ac = A.core;

rca = Ac.r;
rfa = cell(d,1);
rcx = xc.r;
rfx = cell(d,1);
n = cell(d,1);
m = cell(d,1);

L = zeros(d,1);
for i=1:d
    L(i) = xf{i}.d;
    n{i} = Af{i}.n;
    m{i} = Af{i}.m;
    rfa{i} = Af{i}.r;
    rfx{i} = xf{i}.r;
end;

rfy = cell(d,1);
if (isempty(y))
    yc = tt_rand(1, d, 1);
    yf = cell(d,1);
    for i=1:d
        yf{i} = tt_rand(n{i}, L(i), [1; 1*ones(L(i),1)]);
    end;
else
    yc = y.core;
    yf = y.tuck;
end;
rcy = yc.r;
for i=1:d
    rfy{i} = yf{i}.r;
end;

%  keyboard;

phcl = cell(d+1,1); phcl{1} = 1;
phcr = cell(d+1,1); phcr{d+1} = 1;
phfb = cell(d,1);
for i=1:d
    phfb{i} = cell(L(i),1);
    phfb{i}{1} = 1;
end;
phfc = cell(d,1);
phft = cell(d,1);
for i=1:d
    phft{i} = cell(L(i)+1,1);
    phft{i}{L(i)+1} = 1;
end;

rty = zeros(d,1);
rtx = zeros(d,1);
rta = zeros(d,1);


for swp=1:nswp
    dy_max = 0;
    r_max = 0;
    % QRs and phis
    % Factors
    for i=d:-1:1
        cury = yf{i};
        for j=1:L(i)
            cr = cury{j};
            cr = reshape(cr, rfy{i}(j)*n{i}(j), rfy{i}(j+1));
            [cr, rv]=qr(cr, 0);
            if (j<L(i))
                % We are on factor
                cr2 = cury{j+1};
                cr2 = reshape(cr2, rfy{i}(j+1), n{i}(j+1)*rfy{i}(j+2));
                cr2 = rv*cr2;
                rfy{i}(j+1) = size(cr,2);
                cury{j+1} = reshape(cr2, rfy{i}(j+1), n{i}(j+1), rfy{i}(j+2));
            else
                % We are to work with core
                cr2 = yc{i};
                cr2 = permute(cr2, [2, 1, 3]);
                cr2 = reshape(cr2, rfy{i}(j+1), rcy(i)*rcy(i+1));
                cr2 = rv*cr2;
                rfy{i}(j+1) = size(cr,2);
                cr2 = reshape(cr2, rfy{i}(j+1), rcy(i), rcy(i+1));
                yc{i} = permute(cr2, [2, 1, 3]);
            end;
            cr = reshape(cr, rfy{i}(j), n{i}(j), rfy{i}(j+1));
            cury{j} = cr;
            % update bottom phi
            if (j<L(i))
                phfb{i}{j+1} = compute_next_Phi(phfb{i}{j}, cr, Af{i}{j}, xf{i}{j}, 'lr');
            else
                phfc{i} = compute_next_Phi(phfb{i}{j}, cr, Af{i}{j}, xf{i}{j}, 'lr');
            end;
        end;
        yf{i} = cury;
    end;
    % Now for core
    Acr = build_real_core_matrix(Ac, phfc);
    for i=d:-1:2
        rty(i) = rfy{i}(L(i)+1);
        cr = yc{i};
        cr = reshape(cr, rcy(i), rty(i)*rcy(i+1));
        [cr, rv] = qr(cr.', 0);
        cr2 = yc{i-1};
        rty(i-1) = rfy{i-1}(L(i-1)+1);
        cr2 =  reshape(cr2, rcy(i-1)*rty(i-1), rcy(i));
        cr2 = cr2*(rv.');
        rcy(i) = size(cr, 2);
        cr = reshape(cr.', rcy(i), rty(i), rcy(i+1));
        yc{i-1} = reshape(cr2, rcy(i-1), rty(i-1), rcy(i));
        yc{i} = cr;
        % Update right phi
        phcr{i} = compute_next_Phi(phcr{i+1}, cr, Acr{i}, xc{i}, 'rl');
    end;

    % Now, DMRG over factors
    for i=1:d
        rty(i) = rfy{i}(L(i)+1);
        % Convolve core blocks to the last factor blocks
        % And, as usual, it smells like manure
        curx = xf{i}; curxc = xc{i}; rtx(i) = rfx{i}(L(i)+1);
        cura = Af{i}; curac = Ac{i}; rta(i) = rfa{i}(L(i)+1);
        cury = yf{i}; curyc = yc{i};
        curxc = permute(curxc, [2, 1, 3]);
        curxc = reshape(curxc, rtx(i), rcx(i)*rcx(i+1));

        curx = curx*curxc; % Last block is of size m(L), rc1*rc2

        % reshape [last_mode,last_rank] -> [last_mode*last_rank, 1]
        curx = tt_reshape(curx, (curx.n).*[ones(L(i)-1,1); curx.r(L(i)+1)]);
        curyc = permute(curyc, [2, 1, 3]);
        curyc = reshape(curyc, rty(i), rcy(i)*rcy(i+1));        
        cury = cury*curyc; % Last block is of size n(L), rc1*rc2
        % reshape [last_mode,last_rank] -> [last_mode*last_rank, 1]
        cury = tt_reshape(cury, (cury.n).*[ones(L(i)-1,1); cury.r(L(i)+1)]);

        ph1 = phcl{i};
        ph1 = permute(ph1, [1,3,2]);
        ph1 = reshape(ph1, rcy(i)*rcx(i), rca(i));
        ph2 = phcr{i+1};
        ph2 = permute(ph2, [2, 1,3]);
        ph2 = reshape(ph2, rca(i+1), rcy(i+1)*rcx(i+1));
        curac = reshape(curac, rca(i), rta(i)*rca(i+1));
        curac = ph1*curac;
        curac = reshape(curac, rcy(i)*rcx(i)*rta(i), rca(i+1));
        curac = curac*ph2;
        curac = reshape(curac, rcy(i), rcx(i), rta(i), rcy(i+1), rcx(i+1));
        curac = permute(curac, [3, 1, 4, 2, 5]);
        curac = reshape(curac, rta(i), rcy(i)*rcy(i+1)*rcx(i)*rcx(i+1));
        cura = tt_tensor(cura);
        cura = cura*curac; % New last block
        lasta = cura{L(i)};
        lasta = reshape(lasta, rfa{i}(L(i)), n{i}(L(i)), m{i}(L(i)), rcy(i)*rcy(i+1), rcx(i)*rcx(i+1));
        lasta = permute(lasta, [1, 2, 4, 3, 5]);
        cura{L(i)} = reshape(lasta, rfa{i}(L(i)),  n{i}(L(i))*rcy(i)*rcy(i+1)*m{i}(L(i))*rcx(i)*rcx(i+1));
        % New sizes of 1D tt
        curn = n{i}.*[ones(L(i)-1,1); rcy(i)*rcy(i+1)];
        curm = m{i}.*[ones(L(i)-1,1); rcx(i)*rcx(i+1)];
        cura = tt_matrix(cura, curn, curm);

        % Now we are ready to perform mvk =). The last core is nonorth
        for j=L(i):-1:2
            rx1 = rfx{i}(j-1); rx2 = rfx{i}(j);
            ry1 = rfy{i}(j-1); ry2 = rfy{i}(j);
            ra1 = rfa{i}(j-1); ra2 = rfa{i}(j);
            if (j==L(i))
                rx3 = 1; ry3 = 1; ra3 = 1;
            else
                rx3 = rfx{i}(j+1); ra3 = rfa{i}(j+1); ry3 = rfy{i}(j+1);
            end;

            rhs2 = reshape(phft{i}{j+1}, ry3*ra3, rx3);
            x2 = curx{j};
            x2 = reshape(x2, rx2*curm(j), rx3);
            rhs2 = rhs2*(x2.');
            rhs2 = reshape(rhs2, ry3, ra3, rx2, curm(j));
            rhs2 = permute(rhs2, [2, 4, 1, 3]);
            rhs2 = reshape(rhs2, ra3*curm(j), ry3*rx2);
            a2 = cura{j};
            a2 = permute(a2, [2, 1, 4, 3]);
            a2 = reshape(a2, curn(j)*ra2, ra3*curm(j));
            rhs2 = a2*rhs2;
            % We'll need it to compute new phi later
            rhs2 = reshape(rhs2, curn(j), ra2, ry3, rx2);
            rhs2 = permute(rhs2, [1, 3, 2, 4]);

            rhs = reshape(rhs2, curn(j)*ry3*ra2, rx2);
            x1 = curx{j-1};
            x1 = reshape(x1, rx1*curm(j-1), rx2);
            rhs = rhs*(x1.');
            rhs = reshape(rhs, curn(j), ry3, ra2, rx1, curm(j-1));
            rhs = permute(rhs, [5, 3, 4, 1, 2]);
            rhs = reshape(rhs, curm(j-1)*ra2, rx1*curn(j)*ry3);
            a1 = cura{j-1};
            a1 = reshape(a1, ra1*curn(j-1), curm(j-1)*ra2);
            rhs = a1*rhs;
            rhs = reshape(rhs, ra1, curn(j-1), rx1, curn(j)*ry3);
            rhs = permute(rhs, [1, 3, 2, 4]);
            rhs = reshape(rhs, ra1*rx1, curn(j-1)*curn(j)*ry3);
            rhs = reshape(phfb{i}{j-1}, ry1, ra1*rx1)*rhs;
            rhs = reshape(rhs, ry1*curn(j-1), curn(j)*ry3);

            y_prev = cury{j-1};
            y_prev = reshape(y_prev, ry1*curn(j-1), ry2);
            y_prev = y_prev*reshape(cury{j}, ry2, curn(j)*ry3);

            dy = norm(rhs-y_prev, 'fro')/norm(rhs, 'fro');
            dy_max = max(dy_max, dy);

            [u,s,v]=svd(rhs, 'econ');
            s = diag(s);
            nrm = norm(s);
            r = my_chop2(s, eps*nrm/sqrt(L(i))/sqrt(d));
            v = conj(v(:,1:r));
            u = u(:,1:r)*diag(s(1:r));
            % Kick
            v = reort(v, randn(curn(j)*ry3, kickrank));
            radd = size(v,2)-r;
            u = [u, zeros(ry1*curn(j-1), radd)];
            r = r+radd;
            % Stuff back
            cury{j} = reshape(v.', r, curn(j), ry3);
            cury{j-1} = reshape(u, ry1, curn(j-1), r);
            % And update phi
            rhs2 = reshape(rhs2, curn(j)*ry3, ra2*rx2);
            rhs2 = (v')*rhs2;
            phft{i}{j} = reshape(rhs2, r, ra2, rx2);
            rfy{i}(j) = r;
            r_max = max(r_max, r);
            if (verb>1)
                fprintf('=mvrk= swp %d, factor {%d}{%d}, dy: %3.3e, r: %d\n', swp, i, j, dy, r);
            end;
        end;
        % Go back, QR, update phfb, split the tucker core
        for j=1:L(i)
            cr = cury{j};
            if (j<L(i))
                cr = reshape(cr, rfy{i}(j)*n{i}(j), rfy{i}(j+1));
                [cr, rv]=qr(cr, 0);
                cr2 = cury{j+1};
                n2 = size(cr2, 2); ry3 = size(cr2, 3);
                cr2 = reshape(cr2, rfy{i}(j+1), n2*ry3);
                cr2 = rv*cr2;
                rfy{i}(j+1) = size(cr,2);
                cury{j+1} = reshape(cr2, rfy{i}(j+1), n2, ry3);
                cr = reshape(cr, rfy{i}(j), n{i}(j), rfy{i}(j+1));
                cury{j} = cr;
                phfb{i}{j+1} = compute_next_Phi(phfb{i}{j}, cr, Af{i}{j}, xf{i}{j}, 'lr');
            else
                % Split tucker factor
                ry1 = rfy{i}(j); n1 = n{i}(j); n2 = rcy(i)*rcy(i+1);
                cr = reshape(cr, ry1*n1, n2);
                [u,s,v]=svd(cr, 'econ');
                s = diag(s);
                nrm = norm(s);
                r = my_chop2(s, eps*nrm/sqrt(L(i))/sqrt(d));
                u = u(:,1:r);
                v = diag(s(1:r))*(v(:,1:r)');
                % Kick
                u = reort(u, randn(ry1*n1, kickrank));
                radd = size(u,2)-r;
                v = [v; zeros(radd, n2)];
                r = r+radd;
                u = reshape(u, ry1, n1, r);
                cury{j} = u;
                v = reshape(v, r, rcy(i), rcy(i+1));
                yc{i} = permute(v, [2, 1, 3]);
                rfy{i}(L(i)+1) = r;
                r_max = max(r_max, r);
                phfc{i} = compute_next_Phi(phfb{i}{j}, u, Af{i}{j}, xf{i}{j}, 'lr');
                if (verb>1)
                    fprintf('=mvrk= swp %d, tucker_rank(%d), r: %d\n', swp, i, r);
                end;
            end;
        end;
        yf{i} = cury;
        if (i<d)
            % Now, yc{i} is not orthogonal. We can perform a DMRG step for core
            % Acr can be used, except the i-th core, we have to recompute it
            rx1 = rcx(i); rx2 = rcx(i+1); rx3 = rcx(i+2);
            ry1 = rcy(i); ry2 = rcy(i+1); ry3 = rcy(i+2);
            ra1 = rca(i); ra2 = rca(i+1); ra3 = rca(i+2);
            curn1 = rfy{i}(L(i)+1); curn2 = rfy{i+1}(L(i+1)+1);
            curm1 = rfx{i}(L(i)+1); curm2 = rfx{i+1}(L(i+1)+1);
            ph = phfc{i};
            ph = permute(ph, [1, 3,  2]);
            ph = reshape(ph, curn1*curm1, rfa{i}(L(i)+1));
            cura = Ac{i};
            cura = permute(cura, [2, 1, 3]);
            cura = reshape(cura, rfa{i}(L(i)+1), rca(i)*rca(i+1));
            cura = ph*cura;
            cura = reshape(cura, curn1*curm1, rca(i), rca(i+1));
            cura = permute(cura, [2, 1, 3]);
            cura = reshape(cura, rca(i), curn1, curm1, rca(i+1));
            Acr{i} = cura;
            cura = permute(cura, [2, 4, 1, 3]);
            cura = reshape(cura, curn1*ra2, ra1*curm1);

            rhs1 = reshape(phcl{i}, ry1*ra1, rx1);
            rhs1 = rhs1*reshape(xc{i}, rx1, curm1*rx2);
            rhs1 = reshape(rhs1, ry1, ra1, curm1, rx2);
            rhs1 = permute(rhs1, [2, 3, 1, 4]);
            rhs1 = reshape(rhs1, ra1*curm1, ry1*rx2);
            rhs1 = cura*rhs1;
            rhs1 = reshape(rhs1, curn1, ra2, ry1, rx2);
            % This is to be saved
            rhs1 = permute(rhs1, [3, 1, 2, 4]);

            rhs = reshape(rhs1, ry1*curn1*ra2, rx2);
            rhs = rhs*reshape(xc{i+1}, rx2, curm2*rx3);
            rhs = reshape(rhs, ry1*curn1, ra2, curm2, rx3);
            rhs = permute(rhs, [3, 2, 4, 1]);
            rhs = reshape(rhs, curm2*ra2, rx3*ry1*curn1);
            cura = Acr{i+1};
            cura = permute(cura, [2, 4, 3, 1]);
            cura = reshape(cura, curn2*ra3, curm2*ra2);
            rhs = cura*rhs;
            rhs = reshape(rhs, curn2, ra3*rx3, ry1*curn1);
            rhs = permute(rhs, [2, 3, 1]);
            rhs = reshape(rhs, ra3*rx3, ry1*curn1*curn2);
            rhs = reshape(phcr{i+2}, ry3, ra3*rx3)*rhs;
            rhs = reshape(rhs.', ry1*curn1, curn2*ry3);

            y_prev = reshape(yc{i}, ry1*curn1, ry2);
            y_prev = y_prev*reshape(yc{i+1}, ry2, curn2*ry3);

            dy = norm(rhs-y_prev, 'fro')/norm(rhs, 'fro');
            dy_max = max(dy_max, dy);

            [u,s,v]=svd(rhs, 'econ');
            s = diag(s);
            nrm = norm(s);
            r = my_chop2(s, eps*nrm/sqrt(d));
            u = u(:,1:r);
            v = conj(v(:,1:r))*diag(s(1:r));
            % Kick
            u = reort(u, randn(ry1*curn1, kickrank));
            radd = size(u,2)-r;
            v = [v, zeros(curn2*ry3, radd)];
            r = r+radd;
            yc{i} = reshape(u, ry1, curn1, r);
            yc{i+1} = reshape(v.', r, curn2, ry3);
            % Update phi
            rhs1 = reshape(rhs1, ry1*curn1, ra2*rx2);
            phcl{i+1} = reshape((u')*rhs1, r, ra2, rx2);
            rcy(i+1)=r;
            r_max = max(r_max, r);
            if (verb>1)
                fprintf('=mvrk= swp %d, core {%d}, dy: %3.3e, r: %d\n', swp, i, dy, r);
            end;
        end;
    end;

    if (verb>0)
        fprintf('=mvrk= swp %d, dy_max: %3.3e, r_max: %d\n', swp, dy_max, r_max);
    end;
    if (dy_max<eps)
        break;
    end;
end;

y = qtt_tucker;
y.dphys = d;
y.core = yc;
y.tuck = yf;
% y.sz = n;
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

function [Ac]=build_real_core_matrix(Ac, phfc)
% function [Ac]=build_real_core_matrix(Ac, phfc)
% Computes the tt_matrix for core operations from the tucker core of the
% matrix Ac and the factor-core phi matrices
% New matrix sizes are size(phfc{i},1)-by-size(phfc{i},3)

d = Ac.d;
rta = Ac.n;
rca = Ac.r;
rtx = zeros(d,1);
rty = zeros(d,1);

for i=1:d
    ph = phfc{i};
    rtx(i) = size(ph, 1);
    rty(i) = size(ph, 3);
    ph = permute(ph, [1, 3,  2]);
    ph = reshape(ph, rtx(i)*rty(i), rta(i));
    cura = Ac{i};
    cura = permute(cura, [2, 1, 3]);
    cura = reshape(cura, rta(i), rca(i)*rca(i+1));
    cura = ph*cura;
    cura = reshape(cura, rtx(i)*rty(i), rca(i), rca(i+1));
    cura = permute(cura, [2, 1, 3]);
    cura = reshape(cura, rca(i), rtx(i)*rty(i), rca(i+1));
    Ac{i} = cura;
end;

Ac = tt_matrix(Ac, rtx, rty);
end