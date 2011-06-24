function [x]=dmrg_solve2(A, y, x0, eps, tol, rmax, nswp, P, verb)
%  function [x]=dmrg_solve2(A, y, [x0], eps, [tol], [rmax], [nswp], [P],[verb])
% Solves the system P(k,n)*A(n,m)*x(m)=P(k,n)*y(n)
% Default values:
%   x0 = random rank-2
%   tol = eps
%   rmax=1000
%   nswp=15
%   P = I

 

% Inner parameters
max_full_size=2500;
prec_compr=1e-3;
prec_tol=1e-1;
prec_iters=15;
small_verb=true;
use_self_prec=true;
nswp_def=10;
nrestart=40;
gmres_iters=5;
local_prec = 'als';
% local_prec = 'selfprec';
kickrank = 2;


if ((nargin<7)||(isempty(nswp)))
    nswp=nswp_def;
end;
if ((nargin<6)||(isempty(rmax)))
    rmax=1000;
end;
if ((nargin<5)||(isempty(tol)))
    tol=eps;
end;
if ( (nargin<8) || (isempty(verb)) )
  verb=true;
end

input_is_tt_tensor = 0;
if ( isa(A,'tt_matrix') )
  A=core(A);
  input_is_tt_tensor = 1;
  if ((nargin<3)||(isempty(x0)))
      x0=tt_random(tt_size(y), A.tt.d, 2);
  end;
else
    if ((nargin<3)||(isempty(x0)))
        x0=tt_random(tt_size(y), max(size(A)), 2);
    end;
end
if ( nargin >= 8  && isa(P,'tt_matrix') )
  P=core(P);
  input_is_tt_tensor = 1;
end
if ( isa(x0,'tt_tensor') )
  x0=core(x0);
  input_is_tt_tensor = 1;
end
if ( isa(y,'tt_tensor') )
  y=core(y);
  input_is_tt_tensor = 1;
end


d=size(A,1);
if ((nargin<8)||(isempty(P)))
    P = tt_eye(tt_size(y), d);
end;


x=x0;

x{1}=reshape(x{1}, size(x{1},1), 1, size(x{1},2));
y{1}=reshape(y{1}, size(y{1},1), 1, size(y{1},2));
A{1}=reshape(A{1}, size(A{1},1),size(A{1},2), 1, size(A{1},3)); %Bydlocode (@)
P{1}=reshape(P{1}, size(P{1},1),size(P{1},2), 1, size(P{1},3));

phA = cell(d,1);
phy = cell(d,1);

for swp=1:nswp
%     z = x;
    % 1-to-d orthogonalization
    rvx = 1; rnewx=1; phAold=1; phyold=1;

    for i=1:d-1
        cre = x{i};
        n1 = size(cre,1); rx1 = size(cre,2); rx2 = size(cre,3);
        cre = reshape(permute(cre, [2 1 3]), rx1, n1*rx2);
        cre = rvx*cre; % size rnew,n1,rx2
        rx1=rnewx;
        cre = reshape(cre, rx1*n1, rx2);
        [q,rvx]=qr(cre,0); % size rx1*n1,r2new - r2new,rx2
        rnewx = min(rx1*n1, rx2);
        x{i}=permute(reshape(q, rx1, n1, rnewx), [2 1 3]);
        
        % Now, update phi. phA=X' PA X, phY = X' PY 
        a1 = A{i};
        n1=size(a1,1); m1=size(a1,2); ra1=size(a1,3); ra2=size(a1,4);
        p1 = P{i};
        k1=size(p1,1); rp1=size(p1,3); rp2=size(p1,4);
        y1 = y{i};
        ry1=size(y1,2); ry2=size(y1,3);
        x1 = x{i};

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
        phAold=phAold*conj(x1); % size rxm2*rp2*ra2*rxn2 <--- complex conjugate!
        phAold = reshape(phAold, rxm2, rp2, ra2, rxn2);
        phAold = permute(phAold, [4 2 3 1]); % we need rxn2,rp2,ra2,rxm2
        phA{i}=phAold;
        
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
        phyold=phyold*conj(x1); % size rp2*ry2*rxn2 <--- complex conjugate!
        phyold=permute(reshape(phyold, rp2, ry2, rxn2), [3 1 2]);
        phy{i}=phyold;
    end;
    % convolve rv with the last cre
    cre = x{d};
    n1 = size(cre,1); rx1 = size(cre,2); rx2 = size(cre,3);
    cre = reshape(permute(cre, [2 1 3]), rx1, n1*rx2);
    cre = rvx*cre; % size rnew,n1,rx2
    x{d}=permute(reshape(cre, rnewx, n1, rx2), [2 1 3]);
    
    % Now, start the d-to-1 DMRG iteration
    phAold=1; phyold=1;
    for i=d:-1:2
        a2=A{i}; a1=A{i-1}; ra1=size(a1,3); ra2=size(a1,4); ra3=size(a2,4);
        n1 = size(a1,1); m1=size(a1,2); n2=size(a2,1); m2=size(a2,2);
        p2=P{i}; p1=P{i-1}; rp1=size(p1,3); rp2=size(p1,4); rp3=size(p2,4);
        k1 = size(p1,1); k2=size(p2,1);
        
        y1=y{i-1}; y2=y{i}; ry1=size(y1,2); ry2=size(y1,3); ry3=size(y2,3);
        x1=x{i-1}; x2=x{i}; rx1=size(x1,2); rx2=size(x1,3); rx3=size(x2,3);
        
        % Compute RHS: phy{i-2}*P1*y1*y2*P2*phyold
        if (i>2)
            rhs = phy{i-2};
        else
            rhs=1;
        end;
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
        % Now, compress inner rank rp2*ra2 --- what is this ???
        %Modify it by random noise, since sometime MATLAB QR
        %fails
        B2=B2+max(abs(B2(:)))*randn(size(B2))*1e-16;
        [Q,R]=qr(B2,0);

        rnew = min(k2*rxn3*m2*rxm3, rp2*ra2);
        B2 = reshape(Q, k2*rxn3*m2*rxm3, rnew);
        B = B*(R.'); % size rxn1*rxm1*k1*m1*rnew
   
        B = reshape(B, rxn1*k1*rxm1*m1, rnew);
        [U,S,V]=svd(B, 'econ');
        S = diag(S);
        rB = my_chop2(S, 1e-12*norm(S)); % We don't know the cond(B), so let's obtain almost exact compression
        B = U(:,1:rB);
        V = V(:,1:rB)*diag(S(1:rB)); % size rnew*rB        
        B2 = B2*conj(V); % size k2*rxn3*m2*rxm3*rB        
%         rB=rp2*ra2;
        MatVec='bfun2';
        if ((rxn1*k1*k2*rxn3<max_full_size)) % (rB>max(rxn1, rxn3))||
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
        sol_prev = reshape(permute(x{i-1}, [2 1 3]), rxm1*m1, rxm2);
        x2 = reshape(permute(x{i}, [2 1 3]), rxm2, m2*rxm3);
        sol_prev = sol_prev*x2;
        sol_prev = reshape(sol_prev, rxm1*m1*m2*rxm3, 1);

        if (strcmp(MatVec,'full'))

            sol = B \ rhs;
            res=B*sol;

            res_true = norm(res-rhs)/norm(rhs);
        else
            res_prev=norm(bfun2(B,sol_prev,rxm1,m1,m2,rxm3,rxn1,k1,k2,rxn3)-rhs)/norm(rhs);
            [sol_new,flg] = gmres(@(vec)bfun2(B, vec, rxm1,m1,m2,rxm3,rxn1,k1,k2,rxn3), rhs, nrestart, tol, 2, [], [], sol_prev);
            res_new=norm(bfun2(B,sol_new,rxm1,m1,m2,rxm3,rxn1,k1,k2,rxn3)-rhs)/norm(rhs);
            conv_factor=(res_new/res_prev);
            if (res_new*(conv_factor)>eps && use_self_prec) % we need a prec.
%                 keyboard;
                if (strcmp(local_prec, 'selfprec'))
                    iB=tt_minres_selfprec(B, prec_tol, prec_compr, prec_iters, 'right');

                    resid = rhs-bfun2(B,sol_new,rxm1,m1,m2,rxm3,rxn1,k1,k2,rxn3);
                    [dsol,flg] = gmres(@(vec)bfun2(B, bfun2(iB,vec,rxm1,m1,m2,rxm3,rxn1,k1,k2,rxn3),...
                        rxm1,m1,m2,rxm3,rxn1,k1,k2,rxn3), resid, nrestart, tol/res_new, gmres_iters);
                    dsol = bfun2(iB,dsol,rxm1,m1,m2,rxm3,rxn1,k1,k2,rxn3);
                    sol = sol_new+dsol;
                    
                end;
                if (strcmp(local_prec, 'als'))
                    sol = als_solve_rx_2(B, rhs, tol, [], sol_new);
                end;
            else
                [sol,flg] = gmres(@(vec)bfun2(B, vec, rxm1,m1,m2,rxm3,rxn1,k1,k2,rxn3),...
                    rhs, nrestart, tol, gmres_iters, [], [], sol_new);
            end;

            res=bfun2(B,sol,rxm1,m1,m2,rxm3,rxn1,k1,k2,rxn3);
            res_true = norm(res-rhs)/norm(rhs);
        end;
        
        dx = norm(sol-sol_prev,'fro')/norm(sol_prev,'fro');
        if ( verb ) 
        fprintf('==sweep %d, block %d, dx=%3.3e\n', swp, i, dx);
        end
        
        sol=reshape(sol,[rxm1*m1,m2*rxm3]);
        [u,s,v]=svd(sol,'econ');
        s = diag(s);
        flm=norm(s);
        %Truncation block. We have to make it smarter by binary search
        r0=1; r1=min(size(s,1),rmax);
        r=1;
        while ( r ~= r0 || r ~= r1 )
            r=min(floor((r0+r1)/2),rmax);
            er0=norm(s(r+1:numel(s)));
            sol = u(:,1:r)*diag(s(1:r))*(v(:,1:r))';
            sol = reshape(sol, rxm1*m1*m2*rxm3, 1);
            if (strcmp(MatVec,'full'))
                resid = norm(B*sol-rhs)/norm(rhs);
            else
                resid = norm(bfun2(B,sol,rxm1,m1,m2,rxm3,rxn1,k1,k2,rxn3)-rhs)/norm(rhs);
            end;
            if ( verb )
            fprintf('sweep %d, block %d, r=%d, resid=%g, er0=%g, MatVec=%s, rB=%d\n', swp, i, r, resid, er0/flm, MatVec, rB);
            end
            if ((resid<max(res_true*1.2, eps)) ) %Value of the rank is OK
              r1=r;
            else %Is not OK.
              r0=min(r+1,rmax);
            end;
        end
     
        v = conj(v(:,1:r));
        u = u(:,1:r)*diag(s(1:r));
        
        % random kick %This is a production code, sir!
        v = [v, randn(size(v,1), kickrank)];
        u = [u, zeros(size(u,1), kickrank)];
        [v,rv] = qr(v,0);
        r = size(v,2);
        u = u*(rv.');
        
        x{i}=permute(reshape(v, m2, rxm3, r), [1 3 2]);
        x{i-1}=permute(reshape(u, rxm1, m1, r), [2 1 3]);
%         rx2=r;
        rxm2=r; rxn2=r;

        
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
        phAold = phAold*conj(x2); % size rxm2*ra2*rp2*rxn2 <-- cplx conjugate!
        phAold = permute(reshape(phAold, rxm2, ra2, rp2, rxn2), [4 3 2 1]);
        
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
        phyold = phyold*conj(x2); % size rp2*ry2*rxn2 <-- cplx conjugate!
        phyold = permute(reshape(phyold, rp2, ry2, rxn2), [3 1 2]); 
        
    end;
    
end;

x{1}=reshape(x{1}, size(x{1},1), size(x{1},3));

if (input_is_tt_tensor)
  x=tt_tensor(x);
end

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

