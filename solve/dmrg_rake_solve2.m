function [xf, xc]=dmrg_rake_solve2(Af, Ac, yf, yc, tol, varargin)

nswp = 20;
local_format = 'full';
% local_format = 'tt';
max_full_size = 1500;
nrestart = 40;
gmres_iters = 2;
verb = 2;
kickrank = 2;

d = yc.d; % Physical dim.
L = zeros(1,d); % Quantics dims
n = zeros(max(L), d); % Physical mode sizes
for i=1:d
    L(i) = yf{i}.d;
    n(1:L(i), i) = yf{i}.n;
end;

xc = tt_rand(2,d,2);
xf = cell(d,1);
for i=1:d
    xf{i} = tt_rand(n(1:L(i),i), L(i), [1;2*ones(L(i),1)]);
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
phcyl = cell(d+1,1); phcyl{1}=1; % Core, rhs, left
phcyr = cell(d+1,1); phcyr{d+1}=1; % Core, rhs, right
phfyb = cell(d,1); % Factors, rhs, bottom

for swp=1:nswp
    dx_max = 0;
    res_max = 0;
    r_max = 0;
    % bottom-to-top QR and phis
    for i=d:-1:1 % physical dims/core
        phfAb{i} = 1;
        phfyb{i} = 1;
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
            phfAb{i} = compute_next_Phi(phfAb{i}, cr, Af{i}{j}, cr, 'lr');
            phfyb{i} = compute_next_Phi(phfyb{i}, cr, [], yf{i}{j}, 'lr');
        end;
    end;
    
    % DMRG over the core
    % Compute the reduced matrix and rhs;
    rta = zeros(d,1);
    for i=1:d
        rta(i) = rfA(L(i)+1,i);
    end;
    Acr = tt_matrix(Ac, rta, ones(d,1));
    ycr = yc;
    for i=1:d
        a1 = Ac{i};
        a1 = permute(a1, [2, 1, 3]);
        a1 = reshape(a1, rfA(L(i)+1,i), rcA(i)*rcA(i+1));
        ph2 = phfAb{i};
        ph2 = permute(ph2, [1,3,2]);
        ph2 = reshape(ph2, rfx(L(i)+1,i)*rfx(L(i)+1,i), rfA(L(i)+1,i));
        a1 = ph2*a1;
        a1 = reshape(a1, rfx(L(i)+1,i), rfx(L(i)+1,i), rcA(i), rcA(i+1));
        a1 = permute(a1, [3, 1, 2, 4]);
        Acr{i} = a1;
        
        y1 = yc{i};
        y1 = permute(y1, [2, 1, 3]);
        y1 = reshape(y1, rfy(L(i)+1,i), rcy(i)*rcy(i+1));
        ph2 = phfyb{i};
        y1 = ph2*y1;
        y1 = reshape(y1, rfx(L(i)+1,i), rcy(i), rcy(i+1));
        y1 = permute(y1, [2, 1, 3]);
        ycr{i} = y1;        
    end;
    fprintf('=rake_solve========= core processing, sweep %d =======\n', swp);
    xcold = xc;
    xc = dmrg_solve2(Acr, ycr, tol, 'x0', xcold, 'max_full_size', max_full_size, 'nrestart', nrestart, 'gmres_iters', gmres_iters, 'nswp', 1);     
    res = norm(Acr*xcold-ycr)/norm(ycr);
    dx = norm(xcold-xc)/norm(xc);
    dx_max = max(dx_max, dx);
    r_max = max([r_max; rank(xc)]);
    fprintf('=rake_solve========= core res_prev: %3.3e, dx: %3.3e, rmax: %d\n', res, dx, max(rank(xc)));
    res_max = max(res_max, res);
    
    % Horisontal phis
    rcx = xc.r;
    for i=d:-1:2
        % n here is in fact a tucker rank
        cr2 = xc{i};
        cr2 = reshape(cr2, rcx(i), rfx(L(i)+1,i)*rcx(i+1));
        [cr2, rv]=qr(cr2.', 0);
        cr3 = xc{i-1};
        cr3 = reshape(cr3, rcx(i-1)*rfx(L(i-1)+1,i-1), rcx(i));
        cr3 = cr3*(rv.');
        rcx(i) = size(cr2, 2);
        cr2 = cr2.'; % sizes rcx1, rtuck*rcx2
        xc{i} = reshape(cr2, rcx(i), rfx(L(i)+1,i), rcx(i+1));
        xc{i-1} = reshape(cr3, rcx(i-1), rfx(L(i-1)+1,i-1), rcx(i));
        
        % Update right phis
        cr2 = reshape(cr2, rcx(i), rfx(L(i)+1,i), rcx(i+1));
        phcAr{i} = compute_next_Phi(phcAr{i+1}, cr2, Acr{i}, cr2, 'rl');
        % New size: rcx1, rcA1, rcx1
        phcyr{i} = compute_next_Phi(phcyr{i+1}, cr2, [], ycr{i}, 'rl');
        % New size: rcx1, rcy1
    end;
    
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
        a1 = reshape(a1left, rcx(i)*rcx(i)*rta, rcA(i+1))*ph2;
        a1 = reshape(a1, rcx(i), rcx(i), rta, rcx(i+1), rcx(i+1));
        a1 = permute(a1, [1, 4, 2, 5, 3]);
        % we need a1 as well, to check the  residual in tucker block
        % splitting
        a1 = reshape(a1, rcx(i)*rcx(i+1)*rcx(i)*rcx(i+1), rta);
        a2 = Af{i}{L(i)};
        a2 = reshape(a2, rfA(L(i),i)*n(L(i),i)*n(L(i),i), rta);
        a2 = a2*a1.';
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
        y1 = reshape(y1left, rcx(i)*rty, rcy(i+1))*ph2;
        y1 = reshape(y1, rcx(i), rty, rcx(i+1));
        y1 = permute(y1, [1, 3, 2]);
        y1 = reshape(y1, rcx(i)*rcx(i+1), rty);
        y2 = yf{i}{L(i)};
        y2 = reshape(y2, rfy(L(i),i)*n(L(i),i), rty);
        y2 = y2*y1.';
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
        
        fprintf('=rake_solve========= factor {%d} processing, Sweep %d =======\n', i, swp);
        xfrold = xfr;
        xfr = dmrg_solve2(Afr, yfr, tol, 'x0', xfrold, 'max_full_size', max_full_size, 'nrestart', nrestart, 'gmres_iters', gmres_iters, 'nswp', 1);
        res = norm(Afr*xfrold-yfr)/norm(yfr);
        dx = norm(xfrold-xfr)/norm(xfr);
        dx_max = max(dx_max, dx);
        r_max = max([r_max; rank(xfr)]);
        fprintf('=rake_solve========= factor {%d} res_prev: %3.3e, dx: %3.3e, rmax: %d\n', i, res, dx, max(rank(xfr)));
        res_max = max(res_max, res);
        rfx(1:L(i)+1,i) = xfr.r;
        
        % We have to split the tucker block, and compute new phf*b
        phfAb{i}=1;
        phfyb{i}=1;
        for j=1:L(i)
            cr = xfr{j};
            % What is our next core?
            if (j<L(i))
                % we are still on a "tooth"
                cr = reshape(cr, rfx(j,i)*n(j,i), rfx(j+1,i));
                [cr, rv] = qr(cr, 0);
                cr2 = xfr{j+1};
                ncur = size(cr2, 2);
                cr2 = reshape(cr2, rfx(j+1,i), ncur*rfx(j+2,i));
                cr2 = rv*cr2;
                rfx(j+1,i) = size(cr, 2);
                xfr{j} = reshape(cr, rfx(j,i), n(j,i), rfx(j+1,i));
                xfr{j+1} = reshape(cr2, rfx(j+1,i), ncur, rfx(j+2,i));
            else
                % We have to split the tucker core and update the tucker
                % rank
                cr = reshape(cr, rfx(j,i)*n(j,i), rcx(i)*rcx(i+1));
                [u,s,v]=svd(cr, 'econ');
                % Prepare the local matrix and rhs for residue check
                curA = cell(2,1);
                curA{2} = reshape(a1, rcx(i)*rcx(i+1), rcx(i)*rcx(i+1), rta);
                curA{1} = Af{i}{L(i)};
                curA{1} = reshape(curA{1}, rfA(j,i), n(j,i)*n(j,i)*rta);
                ph2 = phfAb{i};
                ph2 = permute(ph2, [1,3,2]);
                ph2 = reshape(ph2, rfx(j,i)*rfx(j,i), rfA(j,i));
                curA{1} = ph2*curA{1};
                curA{1} = reshape(curA{1}, rfx(j,i), rfx(j,i), n(j,i), n(j,i), rta);
                curA{1} = permute(curA{1}, [1, 3, 2, 4, 5]);
                curA{1} = reshape(curA{1}, rfx(j,i)*n(j,i), rfx(j,i)*n(j,i), rta);
                rhs = reshape(yf{i}{j}, rfy(j,i), n(j,i)*rty);
                rhs = phfyb{i}*rhs;
                rhs = reshape(rhs, rfx(j,i)*n(j,i), rty);
                rhs = rhs*(y1.');
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
                    fprintf('=tuckrank= swp %d, fact %d, res: %3.3e, r: %d\n', swp, i, res, r);
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
                rtx = r;
                cr = u;
                xfr{j} = reshape(cr, rfx(j,i), n(j,i), r);                
                v = reshape(v.', r, rcx(i), rcx(i+1));
                v = permute(v, [2, 1, 3]);
                xc{i} = v;
                r_max = max(r_max, r);
            end;
            % Update bottom phis
            cr = reshape(cr, rfx(j,i), n(j,i), rfx(j+1,i));
            phfAb{i} = compute_next_Phi(phfAb{i}, cr, Af{i}{j}, cr, 'lr');
            phfyb{i} = compute_next_Phi(phfyb{i}, cr, [], yf{i}{j}, 'lr');            
        end;
        xf{i} = xfr;
        
        % Now, perform the QR to the next factor, update phc*l as well
        if (i<d)
            cr = xc{i};
            cr = reshape(cr, rcx(i)*rtx, rcx(i+1));
            [cr, rv]=qr(cr, 0);
            cr2 = xc{i+1};
            cr2 = reshape(cr2, rcx(i+1), rfx(L(i+1)+1,i+1)*rcx(i+2));
            cr2 = rv*cr2;
            rcx(i+1) = size(cr, 2);
            xc{i+1} = reshape(cr2, rcx(i+1), rfx(L(i+1)+1,i+1), rcx(i+2));
            xc{i} = reshape(cr, rcx(i), rtx, rcx(i+1));
            
            % We have a1left, y1left of sizes rcx(i)(^2), rtuck, rc*(i+1)
            curph = reshape(a1left, rcx(i)*rcx(i), rta, rcA(i+1));
            curph = permute(curph, [2, 1, 3]);
            curph = reshape(curph, rta, rcx(i)*rcx(i)*rcA(i+1));
            ph2 = phfAb{i};
            ph2 = permute(ph2, [1, 3, 2]);
            ph2 = reshape(ph2, rtx*rtx, rta);
            curph = ph2*curph;
            curph = reshape(curph, rtx, rtx, rcx(i), rcx(i), rcA(i+1));
            curph = permute(curph, [3, 1, 5, 4, 2]);
            curph = reshape(curph, rcx(i)*rtx*rcA(i+1), rcx(i)*rtx);
            curph = curph*cr;
            curph = reshape(curph, rcx(i)*rtx, rcA(i+1)*rcx(i+1));
            curph = (cr')*curph;
            phcAl{i+1} = reshape(curph, rcx(i+1), rcA(i+1), rcx(i+1));
            
            curph = reshape(y1left, rcx(i), rty, rcy(i+1));
            curph = permute(curph, [2, 1, 3]);
            curph = reshape(curph, rty, rcx(i)*rcy(i+1));
            ph2 = phfyb{i};
            curph = ph2*curph;
            curph = reshape(curph, rtx, rcx(i), rcy(i+1));
            curph = permute(curph, [2, 1, 3]);
            curph = reshape(curph, rcx(i)*rtx, rcy(i+1));
            curph = (cr')*curph;
            phcyl{i+1} = curph;            
        end;
    end;
    
    if (verb>0)
        fprintf('\n=rack_solve= swp %d, dx_max: %3.3e, res_max: %3.3e, r_max: %d\n\n', swp, dx_max, res_max, r_max);
    end;
    
    if (res_max<tol)
        break;
    end;
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
