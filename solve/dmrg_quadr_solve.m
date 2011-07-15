function [x,swp]=dmrg_quadr_solve(A, Q, f, x0, eps, tol, rmax, nswp)

maxit = 10;
maxout = 1;
kickrank = 2;
checkrank = 5;
gres = 1;
gmaxit = 20;

if (nargin<6)||(isempty(tol))
    tol=eps;
end;
if (nargin<7)||(isempty(rmax))
    rmax=1000;
end;
if (nargin<8)||(isempty(nswp))
    nswp=1;
end;

n = A.n;
m = A.m;
d = f.d;
ttA = A.tt;
ra = ttA.r;
psa = ttA.ps;
cra = ttA.core;

rf = f.r;
psf = f.ps;
crf = f.core;

ttQ = Q.tt;
rq = ttQ.r;
psq = ttQ.ps;
crq = ttQ.core;

if (nargin<4)||(isempty(x0))
    x0=tt_tensor(tt_ones(d, m));
%     x0 = tt_rand_pos(m,d,1)+tt_tensor(tt_ones(d,m));
end;
checkrank = [1, checkrank*ones(1,d-1), 1];

x = x0;
rx = x.r;
psx = x.ps;
crx = x.core;

phA = cell(d+1,1); phA{1}=1; phA{d+1}=1;
phf = cell(d+1,1); phf{1}=1; phf{d+1}=1;
phQ = cell(d+1,1); phQ{1}=1; phQ{d+1}=1;
for swp=1:nswp
    
    % random checks for residual
    chkA=1; chkQ = 1; chkf = 1;
    chkcr = tt_random(n, d, max(checkrank));
    % Reorthogonalization and Phi
    curpos = psx(d+1);
    for i=d:-1:2
        x1 = crx(psx(i+1)-rx(i)*m(i)*rx(i+1):psx(i+1)-1);
        x1 = reshape(x1, rx(i), m(i)*rx(i+1));
        [x1, rv] = qr(x1.', 0); % m*rx2, rnew - rnew,r1
        rnew = size(x1,2);
        x1 = reshape(x1, m(i), rx(i+1), rnew);
        x1 = permute(x1, [3 1 2]);
        x2 = crx(psx(i-1):psx(i)-1);
        x2 = reshape(x2, rx(i-1)*m(i-1), rx(i));
        x2 = x2*(rv.'); % size r1*m1, rnew
        rx(i)=rnew;
        crx(curpos-rx(i)*m(i)*rx(i+1):curpos-1) = x1(:);
        curpos = curpos-rx(i)*m(i)*rx(i+1);
        psx(i)=curpos;
        crx(curpos-rx(i-1)*m(i-1)*rx(i):curpos-1) = x2(:);        
        
        % Update phA
        a1 = cra(psa(i):psa(i+1)-1);
        a1 = reshape(a1, ra(i), n(i), m(i), ra(i+1));
        phA{i}=reshape(phA{i+1}, rx(i+1)*ra(i+1), rx(i+1));
        curx = reshape(x1, rx(i)*m(i), rx(i+1));
        phA{i}=phA{i}*(curx.'); % size rx2'*ra2, rx1*m1
        phA{i}=reshape(phA{i}, rx(i+1), ra(i+1), rx(i), m(i));
        phA{i}=permute(phA{i}, [1 3 2 4]);
        phA{i}=reshape(phA{i}, rx(i+1)*rx(i), ra(i+1)*m(i));
        a1 = permute(a1, [4 3 2 1]);
        a1 = reshape(a1, ra(i+1)*m(i), n(i)*ra(i));
        phA{i}=phA{i}*a1; % size rx2'*rx1, n1*ra1
        phA{i}=reshape(phA{i}, rx(i+1), rx(i), n(i), ra(i));
        phA{i}=permute(phA{i}, [4 2 3 1]);
        phA{i}=reshape(phA{i}, ra(i)*rx(i), n(i)*rx(i+1));
        curx = reshape(curx, rx(i), m(i)*rx(i+1));
        phA{i}=conj(curx)*(phA{i}.'); % size rx1', ra1*rx1
        phA{i}=reshape(phA{i}, rx(i), ra(i), rx(i));
        
        % Update phQ
        phQ{i}=reshape(phQ{i+1}, rx(i+1)*rq(i+1)*rx(i+1), rx(i+1));
        curx = reshape(curx, rx(i)*m(i), rx(i+1));
        phQ{i}=phQ{i}*(curx.'); % size rx2''*rq2*rx2', rx1*m1
        phQ{i}=reshape(phQ{i}, rx(i+1), rq(i+1), rx(i+1), rx(i), m(i));
        phQ{i}=permute(phQ{i}, [1 2 4 5 3]);
        phQ{i}=reshape(phQ{i}, rx(i+1)*rq(i+1)*rx(i)*m(i), rx(i+1)); %rx2'',rx1,rx2'
        phQ{i}=phQ{i}*(curx.'); % size rx2''*rq2*rx1*m1, rx1'*m1'
        phQ{i}=reshape(phQ{i}, rx(i+1), rq(i+1), rx(i), m(i), rx(i), m(i));
        phQ{i}=permute(phQ{i}, [1 5 3 6 4 2]);
        phQ{i}=reshape(phQ{i}, rx(i+1)*rx(i)*rx(i), m(i)*m(i)*rq(i+1)); %rx2''*rx1'*rx1,m1'*m1*rq2
        q1 = crq(psq(i):psq(i+1)-1);
        q1 = reshape(q1, rq(i)*n(i), m(i)*m(i)*rq(i+1));
        phQ{i}=phQ{i}*(q1.'); % size rx2''*rx1'*rx1, rq1*n1
        phQ{i}=reshape(phQ{i}, rx(i+1),rx(i),rx(i),rq(i),n(i));
        phQ{i}=permute(phQ{i}, [5 1 4 2 3]);
        phQ{i}=reshape(phQ{i}, n(i)*rx(i+1), rq(i)*rx(i)*rx(i));
        curx = reshape(curx, rx(i), m(i)*rx(i+1));
        phQ{i}=conj(curx)*phQ{i}; % size rx1'', rq1*rx1'*rx1
        phQ{i}=reshape(phQ{i}, rx(i), rq(i), rx(i), rx(i));
        
        % Update phf
        f1 = crf(psf(i):psf(i+1)-1);
        f1 = reshape(f1, rf(i)*n(i), rf(i+1));
        phf{i}=phf{i+1}*(f1.'); % size rx2', rf1*n1
        phf{i}=reshape(phf{i}, rx(i+1), rf(i), n(i));
        phf{i}=permute(phf{i}, [3 1 2]);
        phf{i}=reshape(phf{i}, n(i)*rx(i+1), rf(i));
        phf{i}=conj(curx)*phf{i}; % size rx1', rf1
        
        % QR for random check;
        chkcr{i} = permute(chkcr{i}, [1 3 2]);
        chkcr{i} = reshape(chkcr{i}, n(i)*checkrank(i+1), checkrank(i));
        [chkcr{i},rv]=qr(chkcr{i},0);
        chkcr{i-1} = reshape(chkcr{i-1}, n(i-1)*checkrank(i-1), checkrank(i));
        checkrank(i) = size(chkcr{i},2);
        chkcr{i-1}=chkcr{i-1}*(rv.');
        chkcr{i-1} = reshape(chkcr{i-1}, n(i-1), checkrank(i-1), checkrank(i));
        chkcr{i} = reshape(chkcr{i}, n(i), checkrank(i+1), checkrank(i));
        chkcr{i} = permute(chkcr{i}, [1 3 2]);
    end;    
    % Truncate crx.
    truncsize = curpos-rx(1)*m(1)*rx(2)-1;
    crx = crx(truncsize+1:psx(d+1)-1);
    psx(2:d+1)=psx(2:d+1)-truncsize;
    
    res_max = 0;
    % Now, the DMRG scheme
    for i=1:d-1
        % Extract all blocks
        a1 = cra(psa(i):psa(i+1)-1);
        a1 = reshape(a1, ra(i), n(i), m(i), ra(i+1));
        a2 = cra(psa(i+1):psa(i+2)-1);
        a2 = reshape(a2, ra(i+1), n(i+1), m(i+1), ra(i+2));
        f1 = crf(psf(i):psf(i+1)-1);
        f1 = reshape(f1, rf(i), n(i), rf(i+1));
        f2 = crf(psf(i+1):psf(i+2)-1);
        f2 = reshape(f2, rf(i+1), n(i+1), rf(i+2));
        q1 = crq(psq(i):psq(i+1)-1);
        q1 = reshape(q1, rq(i), n(i), m(i), m(i), rq(i+1));
        q2 = crq(psq(i+1):psq(i+2)-1);
        q2 = reshape(q2, rq(i+1), n(i+1), m(i+1), m(i+1), rq(i+2));
        x1 = crx(psx(i):psx(i+1)-1);
        x1 = reshape(x1, rx(i), m(i), rx(i+1));        
        x2 = crx(psx(i+1):psx(i+2)-1);
        x2 = reshape(x2, rx(i+1), m(i+1), rx(i+2));                
        
        % zero-term
        rhs = phf{i}*reshape(f1, rf(i), n(i)*rf(i+1));
        rhs = reshape(rhs, rx(i), n(i), rf(i+1));
        rhs = reshape(rhs, rx(i)*n(i), rf(i+1));
        rhs = rhs*reshape(f2, rf(i+1), n(i+1)*rf(i+2));
        rhs = reshape(rhs, rx(i)*n(i)*n(i+1), rf(i+2));
        rhs = rhs*(phf{i+2}.');
        rhs = reshape(rhs, rx(i)*n(i)*n(i+1)*rx(i+2), 1);
        
        % first-order term (matrix)
        A1 = permute(phA{i}, [1 3 2]);
        A1 = reshape(A1, rx(i)*rx(i), ra(i));
        A1 = A1*reshape(a1, ra(i), n(i)*m(i)*ra(i+1));
        A1 = reshape(A1, rx(i), rx(i), n(i), m(i), ra(i+1));
        A1 = permute(A1, [1 3 2 4 5]);
        A1 = reshape(A1, rx(i)*n(i), rx(i)*m(i), ra(i+1));
        A2 = permute(phA{i+2}, [1 3 2]);
        A2 = reshape(A2, rx(i+2)*rx(i+2), ra(i+2));
        A2 = A2*(reshape(a2, ra(i+1)*n(i+1)*m(i+1), ra(i+2)).');
        A2 = reshape(A2, rx(i+2), rx(i+2), ra(i+1), n(i+1), m(i+1));
        A2 = permute(A2, [4 1 5 2 3]);
        A2 = reshape(A2, n(i+1)*rx(i+2), m(i+1)*rx(i+2), ra(i+1));
        Aloc = cell(2,1);
        Aloc{1}=A1; Aloc{2}=A2;
        Aloc = tt_matrix(Aloc);
        
        % second-order term (3-array)
        Q1 = permute(phQ{i}, [1 3 4 2]);
        Q1 = reshape(Q1, rx(i)*rx(i)*rx(i), rq(i));
        Q1 = Q1*reshape(q1, rq(i), n(i)*m(i)*m(i)*rq(i+1));
        Q1 = reshape(Q1, rx(i), rx(i), rx(i), n(i), m(i), m(i), rq(i+1));
        Q1 = permute(Q1, [1 4 2 5 3 6 7]);
        Q1 = reshape(Q1, rx(i)*n(i), rx(i)*m(i), rx(i)*m(i), rq(i+1));
        Q2 = permute(phQ{i+2}, [1 3 4 2]);
        Q2 = reshape(Q2, rx(i+2)*rx(i+2)*rx(i+2), rq(i+2));
        Q2 = Q2*(reshape(q2, rq(i+1)*n(i+1)*m(i+1)*m(i+1), rq(i+2)).');
        Q2 = reshape(Q2, rx(i+2), rx(i+2), rx(i+2), rq(i+1), n(i+1), m(i+1), m(i+1));
        Q2 = permute(Q2, [4 5 1 6 2 7 3]);
        Q2 = reshape(Q2, rq(i+1), n(i+1)*rx(i+2), m(i+1)*rx(i+2), m(i+1)*rx(i+2));
        Qloc = tt_tensor;
        Qloc.d = 2; Qloc.r=[1; rq(i+1); 1];
        Qloc.n = [(rx(i)*n(i))^3; (n(i+1)*rx(i+2))^3];
        Qloc.ps = cumsum([1; (rx(i)*n(i))^3*rq(i+1); (n(i+1)*rx(i+2))^3*rq(i+1)]);
        Qloc.core = [Q1(:); Q2(:)];
        Qloc = tt_array(Qloc, [rx(i)*n(i); n(i+1)*rx(i+2)]*[1,1,1]);
        
        % previous solution
        sol_prev = reshape(x1, rx(i)*m(i), rx(i+1))*reshape(x2, rx(i+1), m(i+1)*rx(i+2));
        sol_prev = reshape(sol_prev, rx(i)*m(i)*m(i+1)*rx(i+2), 1);
        
%         % Check the previous residual
        res_prev = Aloc*sol_prev + conv3full(Qloc, sol_prev, sol_prev) - rhs;
        res_prev = norm(res_prev)/norm(rhs);
        if (res_prev>res_max)
            res_max = res_prev;
        end;
        
        % Run the newton solver
%         [sol_new,res_new] = cg_2d_solve(rhs, Aloc, Qloc, sol_prev, tol, maxit, maxout);
%         if (res_new>tol)
            [sol_new,res_new] = newton_2d_solve_jac(rhs, Aloc, Qloc, sol_prev, tol, maxit, gres, gmaxit);
%             [sol_new,res_new] = newton_2d_solve(rhs, Aloc, Qloc, sol_new, min(tol,1), 3, gres, gmaxit);
%         end;
        
        sol_new = reshape(sol_new, rx(i)*m(i), m(i+1)*rx(i+2));
        [u,s,v]=svd(sol_new, 'econ');        
%         s = diag(s);
%         dsol = tt_tensor;
%         dsol.d=2;
%         dsol.r = [1;1;1];
%         dsol.n = [rx(i)*m(i); m(i+1)*rx(i+2)];
%         dsol.ps = cumsum([1; rx(i)*m(i); m(i+1)*rx(i+2)]);
%         dsolcr = zeros(rx(i)*m(i)+m(i+1)*rx(i+2), 1);
%         dsol.core=dsolcr;
%         resvec = zeros(numel(sol_new),1);
        % Bin search
        r0 = 1; rM = min(size(s,1),rmax); r = round((r0+rM)/2);
        while (rM-r0>1)
            cur_sol = reshape(u(:,1:r)*s(1:r,1:r)*v(:,1:r)', rx(i)*m(i)*m(i+1)*rx(i+2), 1);
            cur_err = norm(cur_sol(:)-sol_new(:))/norm(sol_new(:));
            cur_res = norm(Aloc*cur_sol + conv3full(Qloc, cur_sol, cur_sol) - rhs)/norm(rhs);
            fprintf('sweep %d, block %d, rank: %d, resid: %3.3e, L2-err: %3.3e\n', swp, i, r, cur_res, cur_err);
%             if ((cur_res<2*res_new)||(cur_res<eps)) && (cur_err<eps)
            if (cur_res<max(2*res_new, eps))
                rM = r-1;
                r = round((r0+rM)/2);
            else
                r0 = r;
                r = round((r0+rM)/2);
            end;
        end;
        % Line search - if the rank is underestimated
        while (r<min(size(s,1), rmax))
            r=r+1;            
%             dsolcr(1:rx(i)*m(i)) = u(:,r);
%             dsolcr(rx(i)*m(i)+1:end) = s(r)*conj(v(:,r));
%             dsol.core=dsolcr;
            cur_sol = reshape(u(:,1:r)*s(1:r,1:r)*v(:,1:r)', rx(i)*m(i)*m(i+1)*rx(i+2), 1);
%             cur_err = norm(s(r+1:end))/norm(s);
            cur_err = norm(cur_sol(:)-sol_new(:))/norm(sol_new(:));
%             resvec = resvec + full(Aloc*dsol + tt_tensor(convten(convten(Qloc,dsol,3), dsol, 2)), numel(sol_new));
%             cur_res = norm(resvec-rhs)/norm(rhs);
            cur_res = norm(Aloc*cur_sol + conv3full(Qloc, cur_sol, cur_sol) - rhs)/norm(rhs);
            fprintf('sweep %d, block %d, rank: %d, resid: %3.3e, L2-err: %3.3e\n', swp, i, r, cur_res, cur_err);
            if ((cur_res<2*res_new)||(cur_res<eps)) && (cur_err<eps)
                break;
            end;
        end;        
%         fprintf('sweep %d, block %d, rank: %d, resid: %3.3e, L2-err: %3.3e\n', swp, i, r, cur_res, cur_err);
        % kickrank
        u = [u(:,1:r), rand(size(u,1), kickrank)];
        v = conj(v(:,1:r))*s(1:r,1:r);
        v = [v, zeros(size(v,1), kickrank)];
        [u,rv]=qr(u,0);
        rx(i+1) = size(u,2);
        v = v*(rv.');
        
        % stuff back. u has already the form r1,n,r2
        v = reshape(v, m(i+1), rx(i+2), rx(i+1));
        v = permute(v, [3 1 2]);
        
        size_x1 = rx(i)*m(i)*rx(i+1); size_x2 = rx(i+1)*m(i+1)*rx(i+2);
        if (size_x1+size_x2>psx(i+2)-psx(i)) 
            % We have to suck and reallocate the storage for x,
            % as our ranks have grown to much
            dsize = size_x1+size_x2 - (psx(i+2)-psx(i));
            crx(psx(i+2)+dsize:psx(d+1)+dsize-1) = crx(psx(i+2):psx(d+1)-1);
            psx(i+2:d+1)=psx(i+2:d+1)+dsize;
        end;
        % Actually stuff u and v
        crx(psx(i):psx(i)+size_x1-1)=u(:);
        psx(i+1)=psx(i)+size_x1;
        crx(psx(i+1):psx(i+1)+size_x2-1)=v(:);
        % Ranks become smaller?
        if (psx(i+2)-(psx(i+1)+size_x2)>0)
            % Shift the core and truncate it
            dsize = psx(i+2)-(psx(i+1)+size_x2);
            crx(psx(i+2)-dsize:psx(d+1)-dsize-1)=crx(psx(i+2):psx(d+1)-1);
            psx(i+2:d+1)=psx(i+2:d+1)-dsize;
            crx = crx(1:psx(d+1)-1);
        end;
        
        % Now, update the Phi and compute the projections for residual
        % checks
        % Current core for random checks
            chkcr{i}=reshape(chkcr{i}, n(i)*checkrank(i), checkrank(i+1));
            [chkcr{i},rv]=qr(chkcr{i},0); % n1*r1, r2new - r2new,r2
            chkcr{i+1}=permute(chkcr{i+1}, [2 1 3]);
            chkcr{i+1}=reshape(chkcr{i+1}, checkrank(i+1), n(i+1)*checkrank(i+2));            
            checkrank(i+1) = size(chkcr{i},2);
            chkcr{i+1}=rv*chkcr{i+1}; % r2new, n2*r3
            chkcr{i+1}=reshape(chkcr{i+1}, checkrank(i+1), n(i+1), checkrank(i+2));
            chkcr{i+1}=permute(chkcr{i+1}, [2 1 3]);
            chkcr{i}=reshape(chkcr{i}, n(i), checkrank(i), checkrank(i+1));
            chkcr{i}=permute(chkcr{i}, [2 1 3]);
            chkcr{i}=reshape(chkcr{i}, checkrank(i)*n(i), checkrank(i+1));
        u = reshape(u, rx(i), m(i)*rx(i+1));
        % Update phA
        a1 = cra(psa(i):psa(i+1)-1);
        a1 = reshape(a1, ra(i), n(i), m(i), ra(i+1));
        phA{i+1}=reshape(phA{i}, rx(i)*ra(i), rx(i));
        phA{i+1}=phA{i+1}*u; % size rx1'*ra1, m1*rx2
            chkA = chkA*u; 
        phA{i+1}=reshape(phA{i+1}, rx(i), ra(i), m(i), rx(i+1));
        phA{i+1}=permute(phA{i+1}, [1 4 2 3]);
        phA{i+1}=reshape(phA{i+1}, rx(i)*rx(i+1), ra(i)*m(i));
            chkA=reshape(chkA, checkrank(i), ra(i), m(i), rx(i+1));
            chkA=permute(chkA, [1 4 2 3]);
            chkA=reshape(chkA, checkrank(i)*rx(i+1), ra(i)*m(i));        
        a1 = permute(a1, [1 3 2 4]);
        a1 = reshape(a1, ra(i)*m(i), n(i)*ra(i+1));
        phA{i+1}=phA{i+1}*a1; % size rx1'*rx2, n1*ra2
            chkA = chkA*a1;
        phA{i+1}=reshape(phA{i+1}, rx(i), rx(i+1), n(i), ra(i+1));
        phA{i+1}=permute(phA{i+1}, [4 2 1 3]);
        phA{i+1}=reshape(phA{i+1}, ra(i+1)*rx(i+1), rx(i)*n(i));
            chkA=reshape(chkA, checkrank(i), rx(i+1), n(i), ra(i+1));
            chkA=permute(chkA, [4 2 1 3]);
            chkA=reshape(chkA, ra(i+1)*rx(i+1), checkrank(i)*n(i));        
        u = reshape(u, rx(i)*m(i), rx(i+1));
        phA{i+1}=(u')*(phA{i+1}.'); % size rx2', ra2*rx2
        phA{i+1}=reshape(phA{i+1}, rx(i+1), ra(i+1), rx(i+1));
            chkA = chkA*chkcr{i};
            chkA = reshape(chkA, ra(i+1), rx(i+1), checkrank(i+1));
            chkA = permute(chkA, [3 1 2]);
            chkA = reshape(chkA, checkrank(i+1)*ra(i+1), rx(i+1));
        
        % Update phQ
        phQ{i+1}=reshape(phQ{i}, rx(i)*rq(i)*rx(i), rx(i));
        u = reshape(u, rx(i), m(i)*rx(i+1));
        phQ{i+1}=phQ{i+1}*u; % size rx1''*rq1*rx1', m1*rx2
            chkQ = chkQ*u;
        phQ{i+1}=reshape(phQ{i+1}, rx(i), rq(i), rx(i), m(i), rx(i+1));
        phQ{i+1}=permute(phQ{i+1}, [1 2 5 4 3]);
        phQ{i+1}=reshape(phQ{i+1}, rx(i)*rq(i)*rx(i+1)*m(i), rx(i)); %rx1'',rx2,rx1'        
            chkQ=reshape(chkQ, checkrank(i), rq(i), rx(i), m(i), rx(i+1));
            chkQ=permute(chkQ, [1 2 5 4 3]);
            chkQ=reshape(chkQ, checkrank(i)*rq(i)*rx(i+1)*m(i), rx(i)); %rx1'',rx2,rx1'        
        phQ{i+1}=phQ{i+1}*u; % size rx1''*rq1*rx2*m1, m1'*rx2'
            chkQ = chkQ*u;
        phQ{i+1}=reshape(phQ{i+1}, rx(i), rq(i), rx(i+1), m(i), m(i), rx(i+1));
        phQ{i+1}=permute(phQ{i+1}, [1 6 3 5 4 2]);
        phQ{i+1}=reshape(phQ{i+1}, rx(i)*rx(i+1)*rx(i+1), m(i)*m(i)*rq(i)); %rx1''*rx2'*rx2,m1'*m1*rq1
            chkQ=reshape(chkQ, checkrank(i), rq(i), rx(i+1), m(i), m(i), rx(i+1));
            chkQ=permute(chkQ, [1 6 3 5 4 2]);
            chkQ=reshape(chkQ, checkrank(i)*rx(i+1)*rx(i+1), m(i)*m(i)*rq(i)); %rx1''*rx2'*rx2,m1'*m1*rq1        
        q1 = crq(psq(i):psq(i+1)-1);
        q1 = reshape(q1, rq(i), n(i), m(i), m(i), rq(i+1));
        q1 = permute(q1, [3 4 1 2 5]);
        q1 = reshape(q1, m(i)*m(i)*rq(i), n(i)*rq(i+1));
        phQ{i+1}=phQ{i+1}*q1; % size rx1''*rx2'*rx2, n1*rq2
            chkQ = chkQ*q1;
        phQ{i+1}=reshape(phQ{i+1}, rx(i),rx(i+1),rx(i+1),n(i),rq(i+1));
        phQ{i+1}=permute(phQ{i+1}, [1 4 5 2 3]);
        phQ{i+1}=reshape(phQ{i+1}, rx(i)*n(i), rq(i+1)*rx(i+1)*rx(i+1));
            chkQ=reshape(chkQ, checkrank(i),rx(i+1),rx(i+1),n(i),rq(i+1));
            chkQ=permute(chkQ, [1 4 5 2 3]);
            chkQ=reshape(chkQ, checkrank(i)*n(i), rq(i+1)*rx(i+1)*rx(i+1));
        u = reshape(u, rx(i)*m(i), rx(i+1));
        phQ{i+1}=(u')*phQ{i+1}; % size rx2'', rq2*rx2'*rx2
        phQ{i+1}=reshape(phQ{i+1}, rx(i+1), rq(i+1), rx(i+1), rx(i+1));
            chkQ = (chkcr{i}.')*chkQ; % size chkrank, rq2*rx2'*rx2
            chkQ = reshape(chkQ, checkrank(i+1)*rq(i+1)*rx(i+1), rx(i+1));
        
        % Update phf
        f1 = crf(psf(i):psf(i+1)-1);
        f1 = reshape(f1, rf(i), n(i)*rf(i+1));
        phf{i+1}=phf{i}*f1; % size rx1', n1*rf2
            chkf = chkf*f1;
        phf{i+1}=reshape(phf{i+1}, rx(i)*n(i), rf(i+1));
            chkf = reshape(chkf, checkrank(i)*n(i), rf(i+1));
        phf{i+1}=(u')*phf{i+1}; % size rx1', rf1        
            chkf = (chkcr{i}.')*chkf;
    end;
    
    % Now, check the residual. We have to convolve Phis with the last cores;
    x1 = crx(psx(d):psx(d+1)-1);
    x1 = reshape(x1, rx(d), n(d));
    a1 = cra(psa(d):psa(d+1)-1);
    a1 = reshape(a1, ra(d), n(d), m(d));
    a1 = permute(a1, [1 3 2]);
    a1 = reshape(a1, ra(d)*m(d), n(d));
    q1 = crq(psq(d):psq(d+1)-1);
    q1 = reshape(q1, rq(d), n(d), m(d), m(d));
    q1 = permute(q1, [3 4 1 2]);
    q1 = reshape(q1, m(d)*m(d)*rq(d), n(d));
    f1 = crf(psf(d):psf(d+1)-1);
    f1 = reshape(f1, rf(d), n(d));
    chkA = chkA*x1; % size chkr*ra1, n(d)
    chkA=reshape(chkA, checkrank(d), ra(d)*m(d));
    chkA = chkA*a1; % size chkr, n(d). chkA is done
    chkQ = chkQ*x1;
    chkQ=reshape(chkQ, checkrank(d), rq(d), rx(d), m(d));
    chkQ=permute(chkQ, [1 2 4 3]);
    chkQ=reshape(chkQ, checkrank(d)*rq(d)*m(d), rx(d)); %rx1'',rx2,rx1'
    chkQ = chkQ*x1;
    chkQ=reshape(chkQ, checkrank(d), rq(d), m(d), m(d));
    chkQ=permute(chkQ, [1 4 3 2]);
    chkQ=reshape(chkQ, checkrank(d), m(d)*m(d)*rq(d)); %rx1''*rx2'*rx2,m1'*m1*rq1
    chkQ = chkQ*q1; % size chkr, n(d). chkQ is done
    chkf = chkf*f1; % size chkr, n(d). chkf is done
    
    res = norm(chkA+chkQ-chkf, 'fro')/norm(chkf, 'fro');    
    
    fprintf('=== sweep %d, check_resid: %3.3e, res_max: %3.3e\n', swp, res, res_max);
    if (res<2*eps)||(res_max<2*eps)
        break;
    end;
%     keyboard;
end;

x.core = crx;
x.r = rx;
x.ps = psx;

end


function [u,err]=newton_2d_solve(rhs, Aloc, Qloc, u0, tol, maxit, gres, gmaxit)
maxout=1;
Qperm2 = permute2(Qloc, [2 3 1]);
Qperm3 = permute2(Qloc, [3 2 1]);
us = zeros(size(rhs,1), maxout);
errs = ones(1,maxout);
for iout=1:maxout
    if (iout==1)
        u=u0;
    else
        u = max(abs(u0))*randn(size(rhs));
    end;
    
    resid = Aloc*u + conv3full(Qloc,u,u) - rhs;
    err_old = norm(resid)/norm(rhs);
    for i=1:maxit
        grad = (Aloc')*resid + conv3full(Qperm2, u, resid) ...
            + conv3full(Qperm3, u, resid);
        
        du = gmres(@(v)HessVec(v, Aloc, Qloc, Qperm2, Qperm3, u, resid), grad, gres, tol, gmaxit);
%         du = bicgstab(@(v)HessVec(v, Aloc, Qloc, Qperm2, Qperm3, u, resid), grad, tol, gmaxit);
        u = u-du;
        resid = Aloc*u + conv3full(Qloc,u,u) - rhs;
        err = norm(resid)/norm(rhs);
        conv = err_old/err;
        err_old = err;
        fprintf('---newton solve: iter [%d,%d], err %3.3e, conv: %3.6f\n',iout,i,err,conv);
        if (err<tol)
            break;
        end;
%         if (conv<1+1e-3)
%             break;
%         end;
    end;
    us(:,iout)=u;
    errs(iout)=err;
    if (err<1e-12)
        break;
    end;
end;
[err,i]=min(errs);
u = us(:,i);
end


function [u,err]=cg_2d_solve(rhs, Aloc, Qloc, u0, tol, maxit, maxout)
if (isempty(tol))
    tol=1e-12;
end;
Qperm2 = permute2(Qloc, [2 3 1]);
Qperm3 = permute2(Qloc, [3 2 1]);
us = zeros(size(rhs,1), maxout);
errs = ones(1,maxout);
for iout=1:maxout
    if (iout==1)
        u=u0;
    else
        u = max(abs(u0))*randn(size(rhs));
    end;
    
    resid = Aloc*u + conv3full(Qloc,u,u) - rhs;
    err_old = norm(resid)/norm(rhs);
    grad_old = zeros(size(rhs));
    dir = zeros(size(rhs));
    for i=1:maxit
        grad = (Aloc')*resid + conv3full(Qperm2, u, resid) ...
            + conv3full(Qperm3, u, resid);
        beta = dot(grad_old,grad_old);
        if (beta~=0)
            beta = dot(grad,grad)/beta;
        end;
        dir = beta*dir+grad;
        
        mix1 = Aloc*dir + conv3full(Qloc, u,dir)+conv3full(Qloc,dir,u);
        mix2 = conv3full(Qloc,dir,dir);
        c0 = dot(resid,resid);
        c1 = 2*dot(mix1,resid);
        c2 = 2*dot(mix2,resid)+dot(mix1,mix1);
        c3 = 2*dot(mix1,mix2);
        c4 = dot(mix2,mix2);
        
        % Frobenius matrix - for roots
        Frob_mtx = zeros(4,4);
        Frob_mtx(2,1)=1;
        Frob_mtx(3,2)=1;
        Frob_mtx(4,3)=1;
        Frob_mtx(1,4)=-c0/c4;
        Frob_mtx(2,4)=-c1/c4;
        Frob_mtx(3,4)=-c2/c4;
        Frob_mtx(4,4)=-c3/c4;
        alpha = eig(Frob_mtx);
        alpha = real(alpha);
        [vv,j]=min(c4*alpha.^4 + c3*alpha.^3 + c2*alpha.^2 + c1*alpha + c0);
        alpha = alpha(j);
%         % A very tiny 1D newton for alpha
%         alpha = -1e-3;
%         for j=1:20
%             agrad = 4*c4*alpha^3 + 3*c3*alpha^2 + 2*c2*alpha + c1;
%             ahess = 12*c4*alpha^2 + 6*c3*alpha + 2*c2;
%             alpha = alpha - agrad/ahess;
%             polyerr=c4*alpha^4 + c3*alpha^3 + c2*alpha^2 + c1*alpha + c0;
%             if (polyerr<1e-13)
%                 break;
%             end;
%         end;
        
        u = u + alpha*dir;
        resid = Aloc*u + conv3full(Qloc,u,u) - rhs;
        err = norm(resid)/norm(rhs);
        conv = err_old/err; 
        err_old = err;
%         fprintf('---cg solve: iter [%d,%d], err %3.3e, conv: %3.6f, poly_err: %3.3e\n',iout,i,err,conv,polyerr);
        grad_old = grad;
        if (err<tol)
            break;
        end;        
        if (conv-1<1e-7)
            break;
        end;
    end;
    fprintf('---cg solve: iter [%d,%d], err %3.3e\n',iout,i,err);
    us(:,iout)=u;
    errs(iout)=err;
    if (err<tol)
        break;
    end;    
end;
[err,i]=min(errs);
u = us(:,i);
end


function [w]=HessVec(v, Aloc, Qloc, Qperm2, Qperm3, u, resid)
w1 = conv3full(Qperm2, v, resid)+conv3full(Qperm3, v, resid);

w2 = Aloc*v + conv3full(Qloc, v, u) + conv3full(Qloc, u, v);

w = (Aloc')*w2 + conv3full(Qperm2, u, w2) ...
    + conv3full(Qperm3, u, w2);
w = w+w1;
end


function [u,err]=newton_2d_solve_jac(rhs, Aloc, Qloc, u0, tol, maxit, gres, gmaxit)
maxout=1;
max_full_size = 200;
if (numel(rhs)<max_full_size)
    Af = full(Aloc);
    Qf = full(Qloc);
end;
us = zeros(size(rhs,1), maxout);
errs = ones(1,maxout);
for iout=1:maxout
    if (iout==1)
        u=u0;
    else
        u = max(abs(u0))*randn(size(rhs));
    end;
    iAloc=[];
    
    resid = Aloc*u + conv3full(Qloc,u,u) - rhs;
    err_old = norm(resid)/norm(rhs);
    for i=1:maxit        
        if (numel(rhs)<max_full_size)
            Bf = reshape(Qf, numel(rhs)*numel(rhs), numel(rhs))*u;
            Bf = reshape(Bf, numel(rhs), numel(rhs));
            Bf2 = reshape(permute(Qf, [1 3 2]), numel(rhs)*numel(rhs), numel(rhs))*u;
            Bf2 = reshape(Bf2, numel(rhs), numel(rhs));
            Bf = Bf + Af + Bf2;
            du = Bf \ resid;
            resgmres=norm(Bf*du-resid)/norm(resid);
        else
            tol_solve = tol*0.1/err_old;
            [du,flg,resgmres_test] = gmres(@(v)JacVec(v, Aloc, Qloc, u), resid, gres, tol_solve, gmaxit);
            u_int = u-du;
            resid_int = Aloc*u_int + conv3full(Qloc,u_int,u_int) - rhs;
            if (norm(resid_int)/norm(rhs)>tol)
%             if (resgmres_test^2>tol_solve)
                if (i==1)
                    iAloc = tt_minres_selfprec(core(Aloc), 1e-1, 1e-3, 10, 'right');
                    iAloc = tt_matrix(iAloc);
                end;
                dres = resid - JacVec(du, Aloc, Qloc, u);
                [ddu,flg,resgmres] = gmres(@(v)JacVec(iAloc*v, Aloc, Qloc, u), dres, gres, tol_solve/resgmres_test, gmaxit);
                du = du + iAloc*ddu;
                resgmres = resgmres*resgmres_test;
            else
                resgmres = resgmres_test;
            end;
        end;
        u = u-du;
        
        resid = Aloc*u + conv3full(Qloc,u,u) - rhs;
        err = norm(resid)/norm(rhs);
        conv = err_old/err;
        err_old = err;
        fprintf('---newton solve: iter [%d,%d], err %3.3e, conv: %3.6f, jac_res: %3.3e\n',iout,i,err,conv, resgmres);
        if (err<tol)
            break;
        end;
        if (conv<1.01)
            break;
        end;
    end;
%     if (~isempty(iAloc))
%         keyboard;
%     end;
    us(:,iout)=u;
    errs(iout)=err;
    if (err<1e-12)
        break;
    end;
end;
[err,i]=min(errs);
u = us(:,i);
end

function [w]=JacVec(v, Aloc, Qloc, u)

w = Aloc*v + conv3full(Qloc,u,v)+conv3full(Qloc,v,u);

end
