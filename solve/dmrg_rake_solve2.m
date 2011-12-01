function [x]=dmrg_rake_solve2(A, y, tol, varargin)

nswp = 20;
local_format = 'full';
% local_format = 'tt';
max_full_size = 1500;
max_full_size2 = 150000;
nrestart = 10;
gmres_iters = 10;
verb = 1;
kickrank = 2;

x = [];

for i=1:2:length(varargin)-1
    switch lower(varargin{i})
        case 'nswp'
            nswp=varargin{i+1};
%         case 'rmax'
%             rmax=lower(varargin{i+1});
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
phcAl = cell(d+1,1); phcAl{1}=1; % Core, matrix, left
phcAr = cell(d+1,1); phcAr{d+1}=1; % Core, matrix, right
phfAb = cell(d,1); % Factors, matrix, bottom
phfAt = cell(d,1); % Factors, matrix, top
phAfc = cell(d,1); % Between the core and the factor
phcyl = cell(d+1,1); phcyl{1}=1; % Core, rhs, left
phcyr = cell(d+1,1); phcyr{d+1}=1; % Core, rhs, right
phfyb = cell(d,1); % Factors, rhs, bottom
phfyt = cell(d,1); % Factors, rhs, top
phyfc = cell(d,1); % Between the core and the factor
for i=1:d
    phfAb{i} = cell(L(i),1);
    phfAb{i}{1}=1;
    phfAt{i} = cell(L(i)+1,1);
    phfAt{i}{L(i)+1}=1;    
    phfyb{i} = cell(L(i),1);
    phfyb{i}{1}=1;
    phfyt{i} = cell(L(i)+1,1);
    phfyt{i}{L(i)+1}=1;        
end;

for swp=1:nswp
    dx_max = 0;
    res_max = 0;
    r_max = 0;
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
                phfAb{i}{j+1} = compute_next_Phi(phfAb{i}{j}, cr, Af{i}{j}, cr, 'lr');
                phfyb{i}{j+1} = compute_next_Phi(phfyb{i}{j}, cr, [], yf{i}{j}, 'lr');
            else
                phAfc{i} = compute_next_Phi(phfAb{i}{j}, cr, Af{i}{j}, cr, 'lr');
                phyfc{i} = compute_next_Phi(phfyb{i}{j}, cr, [], yf{i}{j}, 'lr');                
            end;
        end;
    end;
    
    % QRs and phis over the core
    % Project the system on the factors
    Acr = tt_matrix(Ac, Ac.n, ones(d,1));
    ycr = yc;
    for i=1:d
        Acr{i} = core_matrix(Ac{i}, phAfc{i});
        ycr{i} = core_vector(yc{i}, phyfc{i});
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
        phcAr{i} = compute_next_Phi(phcAr{i+1}, cr, Acr{i}, cr, 'rl');
        phcyr{i} = compute_next_Phi(phcyr{i+1}, cr, [], ycr{i}, 'rl');
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
        rta = rfA(L(i)+1,i); rtx = rfx(L(i)+1,i); rty = rfy(L(i)+1,i);
        a1left = phcAl{i};
        a1left = permute(a1left, [1,3,2]);
        a1left = reshape(a1left, rcx(i)*rcx(i), rcA(i));
        a1left = a1left*reshape(Ac{i}, rcA(i), rta*rcA(i+1));
        % a1left is to use later to compute phcAl
        ph2 = phcAr{i+1};
        ph2 = permute(ph2, [2, 1, 3]);
        ph2 = reshape(ph2, rcA(i+1), rcx(i+1)*rcx(i+1));
        a1top = reshape(a1left, rcx(i)*rcx(i)*rta, rcA(i+1))*ph2;
        a1top = reshape(a1top, rcx(i), rcx(i), rta, rcx(i+1), rcx(i+1));
        a1top = permute(a1top, [1, 4, 2, 5, 3]);
        % we need a1 as well, to check the  residual in tucker block
        % splitting
        a1top = reshape(a1top, rcx(i)*rcx(i+1)*rcx(i)*rcx(i+1), rta);
        a2 = Af{i}{L(i)};
        a2 = reshape(a2, rfA(L(i),i)*n(L(i),i)*n(L(i),i), rta);
        a2 = a2*a1top.';
        a2 = reshape(a2, rfA(L(i),i), n(L(i),i), n(L(i),i), rcx(i)*rcx(i+1), rcx(i)*rcx(i+1));
        a2 = permute(a2, [1, 2, 4, 3, 5]);
        a2 = reshape(a2, rfA(L(i),i), n(L(i),i)*rcx(i)*rcx(i+1), n(L(i),i)*rcx(i)*rcx(i+1), 1);
        Afr = Af{i};
        Afr{L(i)} = a2;
        
        y1left = phcyl{i};
        y1left = y1left*reshape(yc{i}, rcy(i), rty*rcy(i+1));
        % y1left is to use later to compute phcyl
        ph2 = phcyr{i+1};
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
        
        x1 = xc{i};
        x1 = permute(x1, [2, 1, 3]);
        x1 = reshape(x1, rtx, rcx(i)*rcx(i+1));
        x2 = xf{i}{L(i)};
        x2 = reshape(x2, rfx(L(i),i)*n(L(i),i), rtx);
        x2 = x2*x1;
        x2 = reshape(x2, rfx(L(i),i), n(L(i),i)*rcx(i)*rcx(i+1), 1);
        xfr = xf{i};
        xfr{L(i)} = x2;
        
        curn = xfr.n;
        curra = Afr.r;
        curry = yfr.r;
        currx = xfr.r;
        % DMRG over the factor from L(i) to 1
        for j=L(i):-1:2
            a2 = phfAt{i}{j+1};
            a2 = permute(a2, [2, 1, 3]);
            a2 = reshape(a2, curra(j+1), currx(j+1)*currx(j+1));
            a2 = reshape(Afr{j}, curra(j)*curn(j)*curn(j), curra(j+1))*a2;
            a2 = reshape(a2, curra(j), curn(j), curn(j), currx(j+1), currx(j+1));
            a2 = permute(a2, [2, 4, 3, 5, 1]);
            a2 = reshape(a2, curn(j)*currx(j+1), curn(j)*currx(j+1), curra(j));

            a1 = phfAb{i}{j-1};
            a1 = permute(a1, [1, 3, 2]);
            a1 = reshape(a1, currx(j-1)*currx(j-1), curra(j-1));
            a1 = a1*reshape(Afr{j-1}, curra(j-1), curn(j-1)*curn(j-1)*curra(j));
            a1 = reshape(a1, currx(j-1), currx(j-1), curn(j-1), curn(j-1), curra(j));
            a1 = permute(a1, [1, 3, 2, 4, 5]);
            a1 = reshape(a1, currx(j-1)*curn(j-1), currx(j-1)*curn(j-1), curra(j));
            
            y2 = phfyt{i}{j+1}.';            
            y2 = reshape(yfr{j}, curry(j)*curn(j), curry(j+1))*y2;
            y2 = reshape(y2, curry(j), curn(j)*currx(j+1));
            y2 = y2.';
            
            y1 = phfyb{i}{j-1};
            y1 = y1*reshape(yfr{j-1}, curry(j-1), curn(j-1)*curry(j));
            y1 = reshape(y1, currx(j-1)*curn(j-1), curry(j));
            
            x2 = reshape(xfr{j}, currx(j), curn(j)*currx(j+1));
            x2 = x2.';
            x1 = reshape(xfr{j-1}, currx(j-1)*curn(j-1), currx(j));
            
            if (verb>1)
                fprintf('=rake_solve2= swp %d, factor {%d}{%d}, ', swp, i, j);
            end;
            local_format = 'full';
            if (currx(j-1)*curn(j-1)*curn(j)*currx(j+1)>max_full_size2)
                local_format = 'tt';
            end;
            [u,s,v,r,dx_max,res_max]=local_solve(a1, a2, y1, y2, x1, x2, ...
                currx(j-1), curn(j-1), curn(j), currx(j+1), curra(j), ...
                tol/sqrt(L(i))/sqrt(d)/2, res_max, dx_max, ...
                local_format, max_full_size, nrestart, gmres_iters, verb);
            u = u*s;
            % kick
            v = reort(v, randn(curn(j)*currx(j+1), kickrank));
            radd = size(v,2)-r;
            u = [u, zeros(currx(j-1)*curn(j-1), radd)];
            r = r+radd;
            xfr{j} = reshape(v.', r, curn(j), currx(j+1));
            xfr{j-1} = reshape(u, currx(j-1), curn(j-1), r);
            currx(j) = r;
            rfx(j,i)=r;
            r_max = max(r_max, r);
            % new phis
            a2 = permute(a2, [1, 3, 2]);
            a2 = reshape(a2, curn(j)*currx(j+1)*curra(j), curn(j)*currx(j+1));
            a2 = a2*v;
            a2 = reshape(a2, curn(j)*currx(j+1), curra(j)*r);
            a2 = (v')*a2;
            phfAt{i}{j} = reshape(a2, r, curra(j), r);
            phfyt{i}{j} = (v')*y2;
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
                % Prepare the local matrix and rhs for residue check
                curA = cell(2,1);
                curA{2} = reshape(a1top, rcx(i)*rcx(i+1), rcx(i)*rcx(i+1), rta);
                curA{1} = Af{i}{L(i)};
                curA{1} = reshape(curA{1}, rfA(j,i), n(j,i)*n(j,i)*rta);
                ph2 = phfAb{i}{j};
                ph2 = permute(ph2, [1,3,2]);
                ph2 = reshape(ph2, rfx(j,i)*rfx(j,i), rfA(j,i));
                curA{1} = ph2*curA{1};
                curA{1} = reshape(curA{1}, rfx(j,i), rfx(j,i), n(j,i), n(j,i), rta);
                curA{1} = permute(curA{1}, [1, 3, 2, 4, 5]);
                curA{1} = reshape(curA{1}, rfx(j,i)*n(j,i), rfx(j,i)*n(j,i), rta);
                rhs = reshape(yf{i}{j}, rfy(j,i), n(j,i)*rty);
                rhs = phfyb{i}{j}*rhs;
                rhs = reshape(rhs, rfx(j,i)*n(j,i), rty);
                rhs = rhs*(y1top.');
                rhs = reshape(rhs, rfx(j,i)*n(j,i)*rcx(i)*rcx(i+1),1);
                r = 1;
                normy = norm(rhs);
                res_true = norm(bfun2(curA, cr, rfx(j,i), n(j,i), rcx(i), rcx(i+1), rfx(j,i), n(j,i), rcx(i), rcx(i+1))-rhs)/normy;
                while (r<=size(s,1))
                    cursol = u(:,1:r)*s(1:r,1:r)*(v(:,1:r)');
                    res = norm(bfun2(curA, cursol, rfx(j,i), n(j,i), rcx(i), rcx(i+1), rfx(j,i), n(j,i), rcx(i), rcx(i+1))-rhs)/normy;
                    if (res<max(tol/sqrt(L(i)), res_true*2))
                        break;
                    end;
                    r = r+1;
                end;
                if (verb>1)
                    fprintf('=rake_solve2= swp %d, tuckerrank {%d}, res: %3.3e, r: %d\n', swp, i, res, r);
                end;
                u = u(:,1:r);
                v = conj(v(:,1:r));
                s = s(1:r,1:r);
                v = v*s;
                u = reort(u, randn(rfx(j,i)*n(j,i), kickrank));
                radd = size(u,2)-r;
                v = [v, zeros(rcx(i)*rcx(i+1), radd)];
                r = r+radd;
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
                phfAb{i}{j+1} = compute_next_Phi(phfAb{i}{j}, cr, Af{i}{j}, cr, 'lr');
                phfyb{i}{j+1} = compute_next_Phi(phfyb{i}{j}, cr, [], yf{i}{j}, 'lr');            
            else
                phAfc{i} = compute_next_Phi(phfAb{i}{j}, cr, Af{i}{j}, cr, 'lr');
                phyfc{i} = compute_next_Phi(phfyb{i}{j}, cr, [], yf{i}{j}, 'lr');
            end;
        end;
        xf{i} = xfr;
        
        % Now, perform the DMRG over the core to the next factor, update phc*l as well
        if (i<d)
            rtx = rfx(L(i)+1,i); rtx2 = rfx(L(i+1)+1,i+1);
            rx1 = rcx(i); rx2 = rcx(i+1); rx3 = rcx(i+2);
            ra2 = rcA(i+1); ra3 = rcA(i+2); 
            ry2 = rcy(i+1); ry3 = rcy(i+2); 
            % Acr and ycr are okey, except the i-th core
            a1 = reshape(a1left, rx1*rx1, rta, ra2);
            a1 = core_matrix(a1, phAfc{i});
            a1 = reshape(a1, rx1, rx1, rtx, rtx, ra2);
            a1 = permute(a1, [1, 3, 2, 4, 5]);
            a1 = reshape(a1, rx1*rtx, rx1*rtx, ra2);
            
            a2 = phcAr{i+2};
            a2 = permute(a2, [2, 1, 3]);
            a2 = reshape(a2, ra3, rx3*rx3);
            a2 = reshape(Acr{i+1}, ra2*rtx2*rtx2, ra3)*a2;
            a2 = reshape(a2, ra2, rtx2, rtx2, rx3, rx3);
            a2 = permute(a2, [2, 4, 3, 5, 1]);
            a2 = reshape(a2, rtx2*rx3, rtx2*rx3, ra2);
            
            y1 = reshape(y1left, rx1, rty, ry2);
            y1 = core_vector(y1, phyfc{i});
            y1 = reshape(y1, rx1*rtx, ry2);
            
            y2 = phcyr{i+2}.';
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
            [u,s,v,r,dx_max,res_max]=local_solve(a1, a2, y1, y2, x1, x2, ...
                rx1, rtx, rtx2, rx3, ra2, ...
                tol/sqrt(d)/2, res_max, dx_max, ...
                local_format, max_full_size, nrestart, gmres_iters, verb);
            v = v*s;
            % kick
            u = reort(u, randn(rx1*rtx, kickrank));
            radd = size(u,2)-r;
            v = [v, zeros(rtx2*rx3, radd)];
            r = r+radd;
            xc{i} = reshape(u, rx1, rtx, r);
            xc{i+1} = reshape(v.', r, rtx2, rx3);
            rcx(i+1)=r;
            r_max = max(r_max, r);
            % new phis
            a1 = permute(a1, [1, 3, 2]);
            a1 = reshape(a1, rx1*rtx*ra2, rx1*rtx);
            a1 = a1*u;
            a1 = reshape(a1, rx1*rtx, ra2*r);
            a1 = (u')*a1;
            phcAl{i+1} = reshape(a1, r, ra2, r);
            phcyl{i+1} = (u')*y1;            
            
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
        fprintf('\n=rake_solve2= swp %d, dx_max: %3.3e, res_max: %3.3e, r_max: %d\n\n', swp, dx_max, res_max, r_max);
    end;
    
    if (res_max<tol)
        break;
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
ry1 = size(y,1); ry2 = size(y,3);
if (~isempty(A))
  ra1 = size(A,1); ra2 = size(A,4);
else
  ra1 = 1; ra2 = 1;
end

Phi = reshape(Phi_prev, [rx1*ra1, ry1]);
y = reshape(y, [ry1, n*ry2]);
Phi = Phi*y;	% complexity §\mcommentfont$\mathcal{O}(n  r_x r_A r_y^2)$§
Phi = reshape(Phi, [rx1, ra1, n, ry2]);
Phi = permute(Phi, [2, 3, 1, 4]);
if (~isempty(A))
  Phi = reshape(Phi, [ra1*n, rx1*ry2]);
  A = permute(A, [4, 2, 1, 3]);
  A = reshape(A, [ra2*n, ra1*n]);
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
y = B1*x; % size rB*rxn1*k1,m2*rxm3
y = reshape(y, rB, rxn1*k1, m2*rxm3);
y = permute(y, [3 1 2]);
y = reshape(y, m2*rxm3*rB, rxn1*k1);
B2 = reshape(B{2}, k2*rxn3, m2*rxm3*rB);
y = B2*y; % size k2*rxn3,rxn1*k1
y = reshape(y.', rxn1*k1*k2*rxn3, 1);
end


function [u,s,v,r,dx_max,res_max]=local_solve(a1, a2, y1, y2, x1, x2, ...
    rx1, n1, n2, rx3, ra2, ...
    real_tol, res_max, dx_max, ...
    local_format, max_full_size, nrestart, gmres_iters, verb)

if (strcmp(local_format, 'full'))
    sol_prev = x1*(x2.');
    sol_prev = sol_prev(:);
    rhs = y1*(y2.');
    rhs = rhs(:);
    normy = norm(rhs);
    if (rx1*n1*n2*rx3<max_full_size)
        B = reshape(a1, rx1*n1*rx1*n1, ra2);
        B = B*(reshape(a2, n2*rx3*n2*rx3, ra2).');
        B = reshape(B, rx1*n1,  rx1*n1, n2*rx3, n2*rx3);
        B = permute(B, [1, 3, 2, 4]);
        B = reshape(B, rx1*n1*n2*rx3, rx1*n1*n2*rx3);
        sol = B \ rhs;
        res_true = norm(B*sol-rhs)/normy;
        res_prev = norm(B*sol_prev-rhs)/normy;
    else
        B = cell(2,1);
        B{1} = a1;
        B{2} = a2;
        [sol,flg]=gmres(@(v)bfun2(B, v, rx1, n1, n2, rx3, rx1, n1, n2, rx3), rhs, ...
            nrestart, real_tol, gmres_iters, [], [], sol_prev);
        res_true = norm(bfun2(B, sol, rx1, n1, n2, rx3, rx1, n1, n2, rx3)-rhs)/normy;
        res_prev = norm(bfun2(B, sol_prev, rx1, n1, n2, rx3, rx1, n1, n2, rx3)-rhs)/normy;
    end;
    
    dx = norm(sol-sol_prev)/norm(sol);
    dx_max = max(dx_max, dx);
    res_max = max(res_max, res_prev);
    
    sol = reshape(sol, rx1*n1, n2*rx3);
    [u,s,v]=svd(sol, 'econ');
    s = diag(s);
    r = 1;
    while (r<=numel(s))
        cursol = u(:,1:r)*diag(s(1:r))*(v(:,1:r)');
        if (rx1*n1*n2*rx3<max_full_size)
            res = norm(B*cursol(:)-rhs)/normy;
        else
            res = norm(bfun2(B, cursol, rx1, n1, n2, rx3, rx1, n1, n2, rx3)-rhs)/normy;
        end;
        if (res<max(real_tol*2, res_true*2))
            break;
        end;
        r = r+1;
    end;
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
    dsol = tt_gmres(B, drhs, real_tol*2/res_prev, gmres_iters, nrestart, real_tol, real_tol, iB, [], [], [], 1);
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
