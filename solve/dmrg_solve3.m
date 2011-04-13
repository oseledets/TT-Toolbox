function [x]=dmrg_solve3(A, y, x0, eps, tol, rmax, nswp, P)
% function [x]=dmrg_solve3(A, y, [x0], eps, [tol], [rmax], [nswp], [P])
% Solves the system P(k,n)*A(n,m)*x(m)=P(k,n)*y(n)
% Default values:
%   x0 = random rank-2
%   tol = eps
%   rmax=1000
%   nswp=1
%   P = I

% Inner parameters
max_full_size=500;
prec_compr=1e-3;
prec_tol=1e-1;
prec_iters=15;
gmres_iters=20;

d=size(A,1);

if ((nargin<7)||(isempty(nswp)))
    nswp=1;
end;
if ((nargin<6)||(isempty(rmax)))
    rmax=1000;
end;
if ((nargin<5)||(isempty(tol)))
    tol=eps;
end;
if ((nargin<8)||(isempty(P)))
    P = tt_eye(tt_size(y), d);
end;
if ((nargin<3)||(isempty(x0)))
    x0=tt_random(tt_size(y), d, 2);
end;

x=x0;

x{1}=reshape(x{1}, size(x{1},1), 1, size(x{1},2));
y{1}=reshape(y{1}, size(y{1},1), 1, size(y{1},2));
A{1}=reshape(A{1}, size(A{1},1),size(A{1},2), 1, size(A{1},3));
P{1}=reshape(P{1}, size(P{1},1),size(P{1},2), 1, size(P{1},3));

phA = cell(d,1);
phy = cell(d,1);

for swp=1:nswp
%     z = x;
    % 1-to-d orthogonalization
    rvx = 1; rnewx=1; phAold=1; phyold=1;
%     rvz = 1; rnewz=1;
    for i=1:d-1
        core = x{i};
        n1 = size(core,1); rx1 = size(core,2); rx2 = size(core,3);
        core = reshape(permute(core, [2 1 3]), rx1, n1*rx2);
        core = rvx*core; % size rnew,n1,rx2
        rx1=rnewx;
        core = reshape(core, rx1*n1, rx2);
        [q,rvx]=qr(core,0); % size rx1*n1,r2new - r2new,rx2
        rnewx = min(rx1*n1, rx2);
        x{i}=permute(reshape(q, rx1, n1, rnewx), [2 1 3]);
        
%         core = z{i};
%         n1 = size(core,1); rz1 = size(core,2); rz2 = size(core,3);
%         core = reshape(permute(core, [2 1 3]), rz1, n1*rz2);
%         core = rvz*core; % size rnew,n1,rx2
%         rz1=rnewz;
%         core = reshape(core, rz1*n1, rz2);
%         [q,rvz]=qr(core,0); % size rx1*n1,r2new - r2new,rx2
%         rnewz = min(rz1*n1, rz2);
%         z{i}=permute(reshape(q, rz1, n1, rnewz), [2 1 3]);        
        
        % Now, update phi. phA=X^T PA X, phY = X^T PY 
        a1 = A{i};
        n1=size(a1,1); m1=size(a1,2); ra1=size(a1,3); ra2=size(a1,4);
        p1 = P{i};
        k1=size(p1,1); rp1=size(p1,3); rp2=size(p1,4);
        y1 = y{i};
        ry1=size(y1,2); ry2=size(y1,3);
        x1 = x{i};
%         rx1=size(x1,2); rx2=size(x1,3);
%         z1 = z{i};
%         rz1=size(z1,2); rz2=size(z1,3);        
%         phAold = reshape(phAold, rz1*ra1, rx1);
%         x1 = reshape(permute(x1, [2 1 3]), rx1, m1*rx2);
%         phAold=phAold*x1; % size rz1*ra1*m1*rx2
%         phAold=reshape(phAold, rz1, ra1, m1, rx2);
%         phAold=reshape(permute(phAold, [1 4 3 2]), rz1*rx2, m1*ra1);
%         a1 = reshape(permute(a1, [2 3 1 4]), m1*ra1, n1*ra2);
%         phAold=phAold*a1; % size rz1*rx2*n1*ra2
%         phAold=reshape(phAold, rz1, rx2, n1, ra2);
%         phAold=reshape(permute(phAold, [2 4 3 1]), rx2*ra2, n1*rz1);
%         z1 = reshape(z{i}, n1*rz1, rz2);
%         phAold=phAold*z1; % size rx2*ra2*rz2
%         phAold = reshape(phAold, rx2, ra2, rz2);
%         phAold = permute(phAold, [3 2 1]); % we need rz2,ra2,rx2
%         phA{i}=phAold;        
        rxm1=size(x1,2); rxm2=size(x1,3); rxn1=rxm1; rxn2=rxm2;
        phAold = reshape(phAold, rxn1*rp1*ra1, rxm1);
        x1 = reshape(permute(x1, [2 1 3]), rxm1, m1*rxm2);
        phAold=phAold*x1; % size rxn1*rp1*ra1*m1*rxm2
        phAold=reshape(phAold, rxn1, rp1, ra1, m1, rxm2);
        phAold=reshape(permute(phAold, [1 2 5 4 3]), rxn1*rp1*rxm2, m1*ra1);
        a1 = reshape(permute(a1, [2 3 1 4]), m1*ra1, n1*ra2);
        phAold=phAold*a1; % size rxn1*rp1*rxm2*n1*ra2
        phAold=reshape(phAold, rxn1, rp1, rxm2, n1, ra2);
        phAold=reshape(permute(phAold, [1 3 5 4 2]), rxn1*rxm2*ra2, n1*rp1);
        p1 = reshape(permute(p1, [2 3 1 4]), n1*rp1, k1*rp2);
        phAold=phAold*p1; % size rxn1*rxm2*ra2*k1*rp2
        phAold=reshape(phAold, rxn1, rxm2, ra2, k1, rp2);
        phAold=reshape(permute(phAold, [2 5 3 4 1]), rxm2*rp2*ra2, k1*rxn1);
        x1 = reshape(x{i}, k1*rxn1, rxn2);
        phAold=phAold*x1; % size rxm2*rp2*ra2*rxn2
        phAold = reshape(phAold, rxm2, rp2, ra2, rxn2);
        phAold = permute(phAold, [4 2 3 1]); % we need rxn2,rp2,ra2,rxm2
        phA{i}=phAold;
        
%         phyold = reshape(phyold, rz1, ry1);
%         y1 = reshape(permute(y1, [2 1 3]), ry1, n1*ry2);
%         phyold=phyold*y1; % size rz1*n1*ry2
%         phyold = reshape(phyold, rz1, n1, ry2);
%         phyold=reshape(permute(phyold, [3 2 1]), ry2, n1*rz1);
%         z1=reshape(z{i}, n1*rz1, rz2);
%         phyold=phyold*z1; % size ry2*rz2
%         phyold=permute(reshape(phyold, ry2, rz2), [2 1]);
%         phy{i}=phyold;
        phyold = reshape(phyold, rxn1*rp1, ry1);
        y1 = reshape(permute(y1, [2 1 3]), ry1, n1*ry2);
        phyold=phyold*y1; % size rxn1*rp1*n1*ry2
        phyold = reshape(phyold, rxn1, rp1, n1, ry2);
        phyold=reshape(permute(phyold, [1 4 3 2]), rxn1*ry2, n1*rp1);
        p1=reshape(permute(P{i}, [2 3 1 4]), n1*rp1, k1*rp2);
        phyold=phyold*p1; % size rxn1*ry2*k1*rp2
        phyold=reshape(phyold, rxn1, ry2, k1, rp2);
        phyold=reshape(permute(phyold, [4 2 3 1]), rp2*ry2, k1*rxn1);
        x1=reshape(x{i}, k1*rxn1, rxn2);
        phyold=phyold*x1; % size rp2*ry2*rxn2
        phyold=permute(reshape(phyold, rp2, ry2, rxn2), [3 1 2]);
        phy{i}=phyold;
    end;
    % convolve rv with the last core
    core = x{d};
    n1 = size(core,1); rx1 = size(core,2); rx2 = size(core,3);
    core = reshape(permute(core, [2 1 3]), rx1, n1*rx2);
    core = rvx*core; % size rnew,n1,rx2
    x{d}=permute(reshape(core, rnewx, n1, rx2), [2 1 3]);
%     core = z{d};
%     n1 = size(core,1); rz1 = size(core,2); rz2 = size(core,3);
%     core = reshape(permute(core, [2 1 3]), rz1, n1*rz2);
%     core = rvz*core; % size rnew,n1,rx2
%     z{d}=permute(reshape(core, rnewz, n1, rz2), [2 1 3]);    
    
    % Now, start the d-to-1 DMRG iteration
    phAold=1; phyold=1;
    for i=d:-1:2
        a2=A{i}; a1=A{i-1}; ra1=size(a1,3); ra2=size(a1,4); ra3=size(a2,4);
        n1 = size(a1,1); m1=size(a1,2); n2=size(a2,1); m2=size(a2,2);
        p2=P{i}; p1=P{i-1}; rp1=size(p1,3); rp2=size(p1,4); rp3=size(p2,4);
        k1 = size(p1,1); k2=size(p2,1);
        
        y1=y{i-1}; y2=y{i}; ry1=size(y1,2); ry2=size(y1,3); ry3=size(y2,3);
        x1=x{i-1}; x2=x{i}; rx1=size(x1,2); rx2=size(x1,3); rx3=size(x2,3);
%         z1=z{i-1}; z2=z{i}; rz1=size(z1,2); rz3=size(z2,3);
        
        % Compute RHS: phy{i-2}*P1*y1*y2*P2*phyold
        if (i>2)
            rhs = phy{i-2};
        else
            rhs=1;
        end;
%         rhs = reshape(rhs, rz1, ry1);
%         y1 = reshape(permute(y1, [2 1 3]), ry1, n1*ry2);
%         rhs = rhs*y1; % size rz1*n1*ry2
%         rhs = reshape(rhs, rz1*n1, ry2);
%         y2=reshape(permute(y2, [2 1 3]), ry2, n2*ry3);
%         rhs = rhs*y2; % size rz1*n1*n2*ry3
%         rhs = reshape(rhs, rz1*n1*n2, ry3);
%         rhs = rhs*(phyold.'); % size rz1*n1*n2*rz3
%         rhs = reshape(rhs, rz1*n1*n2*rz3, 1);        
        rhs = reshape(rhs, rx1*rp1, ry1);
        y1 = reshape(permute(y1, [2 1 3]), ry1, n1*ry2);
        rhs = rhs*y1; % size rx1*rp1*n1*ry2
        rhs = reshape(rhs, rx1, rp1, n1, ry2);
        rhs=reshape(permute(rhs, [1 4 3 2]), rx1*ry2, n1*rp1);
        p1 = reshape(permute(p1, [2 3 1 4]), n1*rp1, k1*rp2);
        rhs=rhs*p1; % size rx1*ry2*k1*rp2
        rhs=reshape(rhs, rx1, ry2, k1, rp2);
        rhs=reshape(permute(rhs, [1 4 3 2]), rx1*rp2*k1, ry2);
        y2=reshape(permute(y2, [2 1 3]), ry2, n2*ry3);
        rhs = rhs*y2; % size rx1*rp2*k1*n2*ry3
        rhs = reshape(rhs, rx1, rp2, k1, n2, ry3);
        rhs = reshape(permute(rhs, [1 3 5 4 2]), rx1*k1*ry3, n2*rp2);
        p2 = reshape(permute(p2, [2 3 1 4]), n2*rp2, k2*rp3);
        rhs = rhs*p2; % size rx1*k1*ry3*k2*rp3
        rhs = reshape(rhs, rx1, k1, ry3, k2, rp3);
        rhs = reshape(permute(rhs, [1 2 4 5 3]), rx1*k1*k2, rp3*ry3);
        phyold2 = reshape(phyold, rx3, rp3*ry3);
        rhs = rhs*(phyold2.'); % size rx1*k1*k2*rx3
        rhs = reshape(rhs, rx1*k1*k2*rx3, 1);
        
        % Form matrix: phA{i-2}*P1*A1*A2*P2*phAold
%         if (i>2)
%             B = phA{i-2};
%         else
%             B=1;
%         end;
%         B = reshape(B, rz1, ra1, rx1);
%         B = reshape(permute(B, [1 3 2]), rz1*rx1, ra1);
%         a1 = reshape(permute(A{i-1}, [3 2 1 4]), ra1, m1*n1*ra2);
%         B = B*a1; % size rz1*rx1*m1*n1*ra2
%         B = reshape(B, rz1,rx1,m1,n1,ra2);
%         B = permute(B, [1 4 2 3 5]);
%         B = reshape(B, rz1*n1*rx1*m1, ra2);
% %         This is the first term of tensor-structured matrix B \otimes B2
% %         Now, the second
%         B2 = permute(phAold, [1 3 2]);
%         B2 = reshape(B2, rz3*rx3, ra3);
%         a2 = reshape(A{i}, n2*m2*ra2, ra3);
%         B2 = B2*(a2.'); % size rz3*rx3*n2*m2*ra2
%         B2 = reshape(B2, rz3, rx3, n2, m2, ra2);
%         B2 = permute(B2, [3 1 4 2 5]);
%         B2 = reshape(B2, n2*rz3*m2*rx3, ra2);      
%         MatVec='bfun2';
% %         if (ra2*max(rz1,rz3)>rx1*rx3)
% %             MatVec='full';
% %             B = B*(B2.'); % size rz1*n1*rx1*m1*n2*rz3*m2*rx3
% %             B = reshape(B, rz1,n1,rx1,m1,n2,rz3,m2,rx3);
% %             B = permute(B, [1 2 5 6 3 4 7 8]);
% %             B = reshape(B, rz1*n1*n2*rz3, rx1*m1*m2*rx3);
% %         else
%             B = reshape(B, rz1*n1, rx1*m1, ra2);
%             B2 = reshape(B2, n2*rz3, m2*rx3, ra2);
% %         end;        
        rxn1=rx1; rxn3=rx3;
        rxm1=rx1; rxm2=rx2; rxm3=rx3;
        if (i>2)
            B = phA{i-2};
        else
            B=1;
        end;
        B = reshape(B, rxn1, rp1, ra1, rxm1);
        B = reshape(permute(B, [1 4 2 3]), rxn1*rxm1*rp1, ra1);
        a1 = reshape(permute(A{i-1}, [3 2 1 4]), ra1, m1*n1*ra2);
        B = B*a1; % size rxn1*rxm1*rp1*m1*n1*ra2
        B = reshape(B, rxn1,rxm1,rp1,m1,n1,ra2);
        B = permute(B, [1 2 4 6 3 5]);
        B = reshape(B, rxn1*rxm1*m1*ra2, rp1*n1);
        p1 = reshape(permute(P{i-1}, [3 2 1 4]), rp1*n1, k1*rp2);
        B = B*p1; % size rxn1*rxm1*m1*ra2*k1*rp2
        B = reshape(B, rxn1,rxm1,m1,ra2,k1,rp2);
        B = permute(B, [1 5 2 3 6 4]);
        B = reshape(B, rxn1*k1*rxm1*m1, rp2*ra2);
%         This is the first term of tensor-structured matrix B \otimes B2
%         Now, the second
        B2 = permute(phAold, [1 2 4 3]);
        B2 = reshape(B2, rxn3*rp3*rxm3, ra3);
        a2 = reshape(A{i}, n2*m2*ra2, ra3);
        B2 = B2*(a2.'); % size rxn3*rp3*rxm3*n2*m2*ra2
        B2 = reshape(B2, rxn3, rp3, rxm3, n2, m2, ra2);
        B2 = permute(B2, [1 3 5 6 4 2]);
        B2 = reshape(B2, rxn3*rxm3*m2*ra2, n2*rp3);
        p2 = reshape(permute(P{i}, [2 4 1 3]), n2*rp3, k2*rp2);
        B2 = B2*p2; % size rxn3*rxm3*m2*ra2*k2*rp2
        B2 = reshape(B2, rxn3, rxm3, m2, ra2, k2, rp2);
        B2 = permute(B2, [5 1 3 2 6 4]);
        B2 = reshape(B2, k2*rxn3*m2*rxm3, rp2*ra2);
        % Now, compress inner rank rp2*ra2
        [Q,R]=qr(B2,0);
%        if (numel(find(isnan(Q)))>0)
%          keyboard;
%        end

        rnew = min(k2*rxn3*m2*rxm3, rp2*ra2);
        B2 = reshape(Q, k2*rxn3*m2*rxm3, rnew);
        B = B*(R.'); % size rxn1*rxm1*k1*m1*rnew
   
        B = reshape(B, rxn1*k1*rxm1*m1, rnew);
        [U,S,V]=svd(B, 'econ');
        S = diag(S);
        rB = my_chop2(S, eps*norm(S));
        B = U(:,1:rB);
        V = V(:,1:rB)*diag(S(1:rB)); % size rnew*rB        
        B2 = B2*V; % size k2*rxn3*m2*rxm3*rB        
%         rB=rp2*ra2;
        MatVec='bfun2';
        if ((rB>max(rxn1, rxn3))||(rxn1*k1*k2*rxn3<max_full_size))
            MatVec='full';
            B = B*(B2.'); % size rxn1*k1*rxm1*m1*k2*rxn3*m2*rxm3
            B = reshape(B, rxn1,k1,rxm1,m1,k2,rxn3,m2,rxm3);
            B = permute(B, [1 2 5 6 3 4 7 8]);
            B = reshape(B, rxn1*k1*k2*rxn3, rxm1*m1*m2*rxm3);
        else
            B1 = reshape(B, rxn1*k1, rxm1*m1, rB);
            B2 = reshape(B2, k2*rxn3, m2*rxm3, rB);
            B=cell(2,1);
            B{1}=B1;
            B{2}=B2;         
        end;
        
        % Form previous solution
%         sol_prev = reshape(permute(x{i-1}, [2 1 3]), rx1*m1, rx2);
%         x2 = reshape(permute(x{i}, [2 1 3]), rx2, m2*rx3);
%         sol_prev = sol_prev*x2;
%         sol_prev = reshape(sol_prev, rx1*m1*m2*rx3, 1);        
        sol_prev = reshape(permute(x{i-1}, [2 1 3]), rxm1*m1, rxm2);
        x2 = reshape(permute(x{i}, [2 1 3]), rxm2, m2*rxm3);
        sol_prev = sol_prev*x2;
        sol_prev = reshape(sol_prev, rxm1*m1*m2*rxm3, 1);

%         if (strcmp(MatVec,'full'))
%             sol = pcg(B'*B, B'*rhs, tol, 1000, [], [], sol_prev);
%             res_true = norm(B*sol-rhs)/norm(rhs);
%         else
%             rhs2 = bfun3(B, B2, rhs, rx1, m1, m2, rx3, rz1, n1, n2, rz3, ra2, 'BT');
%             sol = pcg(@(vec)bfun3(B, B2, vec, rx1, m1, m2, rx3, rz1, n1, n2, rz3, ra2, 'BTB'), rhs2, tol, 10000, [], [], sol_prev);
%             res_true = norm(bfun3(B, B2, sol, rx1, m1, m2, rx3, rz1, n1, n2, rz3, ra2, 'B')-rhs)/norm(rhs);
%         end;
        if (strcmp(MatVec,'full'))
%             grid_data=amg_grids_setup(B);
%             smooth_params=amg_smoother_params(grid_data);
%             smooth_data=amg_smoother_setup(grid_data, smooth_params);
%             sol_prev2=B*sol_prev;
%             sol = gmres(@(v)(B*amg_ncycles(v,grid_data,smooth_data,5)), rhs, 50, tol, 10, [], [], sol_prev2);
%             sol=amg_ncycles(sol,grid_data,smooth_data,5);
%             sol = gmres(B, rhs, 50, tol, 10, [], [], sol_prev);
            sol = B \ rhs;
            res=B*sol;

            res_true = norm(res-rhs)/norm(rhs);
        else
            %res_prev=norm(bfun2(B,sol_prev,rxm1,m1,m2,rxm3,rxn1,k1,k2,rxn3)-rhs)/norm(rhs);
            %sol_new = gmres(@(vec)bfun2(B, vec, rxm1,m1,m2,rxm3,rxn1,k1,k2,rxn3), rhs, 50, tol, 5, [], [], sol_prev);
            sol = als_solve3(B,rhs,tol,sol_prev,100);
            res=bfun2(B,sol,rxm1,m1,m2,rxm3,rxn1,k1,k2,rxn3);
%             toc;
            res_true = norm(res-rhs)/norm(rhs);
          %  keyboard;
        end;
        
        sol=reshape(sol,[rxm1*m1,m2*rxm3]);
        [u,s,v]=svd(sol,'econ');
        s = diag(s);
        flm=norm(s);
        for r=1:min(size(s,1), rmax)
            er0=norm(s(r+1:numel(s)));
            sol = u(:,1:r)*diag(s(1:r))*(v(:,1:r))';
            sol = reshape(sol, rxm1*m1*m2*rxm3, 1);
            if (strcmp(MatVec,'full'))
                resid = norm(B*sol-rhs)/norm(rhs);
            else
                resid = norm(bfun2(B,sol,rxm1,m1,m2,rxm3,rxn1,k1,k2,rxn3)-rhs)/norm(rhs);
%                 resid = norm(bfun3(B, B2, sol, rx1, m1, m2, rx3, rz1, n1, n2, rz3, ra2, 'B')-rhs)/norm(rhs);
            end;
            fprintf('sweep %d, block %d, r=%d, resid=%g, er0=%g, MatVec=%s, rB=%d\n', swp, i, r, resid, er0/flm, MatVec, rB);
            if ((resid<max(res_true*1.2, eps)) )
                break;
            end;
        end;
        v = v(:,1:r);
        u = u(:,1:r)*diag(s(1:r));
        x{i}=permute(reshape(v, m2, rxm3, r), [1 3 2]);
        x{i-1}=permute(reshape(u, rxm1, m1, r), [2 1 3]);
%         rx2=r;
        rxm2=r; rxn2=r;
%         % Orthogonalize Z
%         core = z{i};
%         n2 = size(core,1); rz2 = size(core,2); rz3 = size(core,3);
%         core = reshape(permute(core, [1 3 2]), n2*rz3, rz2);
%         [q,rvz]=qr(core,0); % size n2,rz3,rz2new - rz2new,rz2
%         rnewz = min(n2*rz3, rz2);
%         z{i}=permute(reshape(q, n2, rz3, rnewz), [1 3 2]);
%         core = z{i-1};
%         n1 = size(core,1); rz1 = size(core,2); rz2 = size(core,3); 
%         core = reshape(core, n1*rz1, rz2);
%         core = core*(rvz.'); % size n1*rz1*rnewz
%         z{i-1}=reshape(core, n1, rz1, rnewz);
%         rz2 = rnewz;
        
        % update phi
%         phAold = reshape(phAold, rz3*ra3, rx3);
%         x2 = reshape(x{i}, m2*rx2, rx3);
%         phAold = phAold*(x2.'); % size rz3*ra3*m2*rx2
%         phAold = reshape(phAold, rz3, ra3, m2, rx2);
%         phAold = permute(phAold, [1 4 3 2]);
%         phAold = reshape(phAold, rz3*rx2, m2*ra3);
%         a2 = reshape(permute(A{i}, [2 4 1 3]), m2*ra3, n2*ra2);
%         phAold=phAold*a2; % size rz3*rx2*n2*ra2
%         phAold=reshape(phAold, rz3, rx2, n2, ra2);
%         phAold = permute(phAold, [2 4 1 3]);
%         phAold = reshape(phAold, rx2*ra2, rz3*n2);        
%         z2 = reshape(permute(z{i}, [3 1 2]), rz3*n2, rz2);
%         phAold = phAold*z2; % size rx2*ra2*rz2
%         phAold = permute(reshape(phAold, rx2, ra2, rz2), [3 2 1]);      
        phAold = reshape(phAold, rxn3*rp3*ra3, rxm3);
        x2 = reshape(x{i}, m2*rxm2, rxm3);
        phAold = phAold*(x2.'); % size rxn3*rp3*ra3*m2*rxm2
        phAold = reshape(phAold, rxn3, rp3, ra3, m2, rxm2);
        phAold = permute(phAold, [1 2 5 4 3]);
        phAold = reshape(phAold, rxn3*rp3*rxm2, m2*ra3);
        a2 = reshape(permute(A{i}, [2 4 1 3]), m2*ra3, n2*ra2);
        phAold=phAold*a2; % size rxn3*rp3*rxm2*n2*ra2
        phAold=reshape(phAold, rxn3, rp3, rxm2, n2, ra2);
        phAold = permute(phAold, [1 3 5 4 2]);
        phAold = reshape(phAold, rxn3*rxm2*ra2, n2*rp3);
        p2 = reshape(permute(P{i}, [2 4 1 3]), n2*rp3, k2*rp2);
        phAold=phAold*p2; % size rxn3*rxm2*ra2*k2*rp2
        phAold=reshape(phAold, rxn3,rxm2,ra2,k2,rp2);
        phAold = permute(phAold, [2 3 5 1 4]);
        phAold = reshape(phAold, rxm2*ra2*rp2, rxn3*k2);        
        x2 = reshape(permute(x{i}, [3 1 2]), rxn3*k2, rxn2);
        phAold = phAold*x2; % size rxm2*ra2*rp2*rxn2
        phAold = permute(reshape(phAold, rxm2, ra2, rp2, rxn2), [4 3 2 1]);
        
%         phyold = reshape(phyold, rz3, ry3);
%         y2 = reshape(y{i}, n2*ry2, ry3);
%         phyold = phyold*(y2.'); % size rz3*n2*ry2
%         phyold = reshape(phyold, rz3, n2, ry2);
%         phyold = permute(phyold, [3 1 2]);
%         phyold = reshape(phyold, ry2, rz3*n2);       
%         z2 = reshape(permute(z{i}, [3 1 2]), rz3*n2, rz2);
%         phyold = phyold*z2; % size ry2*rz2
%         phyold = permute(reshape(phyold, ry2, rz2), [2 1]); 
        phyold = reshape(phyold, rxn3*rp3, ry3);
        y2 = reshape(y{i}, n2*ry2, ry3);
        phyold = phyold*(y2.'); % size rxn3*rp3*n2*ry2
        phyold = reshape(phyold, rxn3, rp3, n2, ry2);
        phyold = permute(phyold, [1 4 3 2]);
        phyold = reshape(phyold, rxn3*ry2, n2*rp3);
        p2 = reshape(permute(P{i}, [2 4 1 3]), n2*rp3, k2*rp2);
        phyold = phyold*p2; % size rxn3*ry2*k2*rp2
        phyold = reshape(phyold, rxn3, ry2, k2, rp2);
        phyold = permute(phyold, [4 2 1 3]);
        phyold = reshape(phyold, rp2*ry2, rxn3*k2);       
        x2 = reshape(permute(x{i}, [3 1 2]), rxn3*k2, rxn2);
        phyold = phyold*x2; % size rp2*ry2*rxn2
        phyold = permute(reshape(phyold, rp2, ry2, rxn2), [3 1 2]); 
        
%         mesh(reshape(tc_to_full(x), 1024*[1,1]));
%         pause(0.5);
    end;
    
end;

x{1}=reshape(x{1}, size(x{1},1), size(x{1},3));

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
