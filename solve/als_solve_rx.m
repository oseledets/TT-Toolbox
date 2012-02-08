function [x]=als_solve_rx(mat, rhs, tol, drx, nswp, addswp)
%Computes an approximate low-rank solution for 2D case
%   [x]=ALS_SOLVE_RX(MAT, RHS, [TOL], [RX], [NSWP])
%   Finds a solution to 2D TTM matrix MAT using the ALS to a 2D TT tensor
%   with rank rx, but the RHS and X are represented as full vectors.
%   TOL is the tolerance for ||x_{i+1}-x_i||/||x_i||,
%   DRX is the random kick rank,
%   NSWP - number of ALS sweeps.
%   default values:
%   tol: 1e-12
%   rx: 1
%   nswp: 10
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


a1 = mat{1}; a2 = mat{2};
n1 = size(a1,1); m1 = size(a1,2);
n2 = size(a2,1); m2 = size(a2,2);
ra = size(a1,3);

% tol2 = 1e-3;

rhs = reshape(rhs, n1, n2);

if (nargin<3)||(isempty(tol))
    tol=1e-12;
end;
if (nargin<4)||(isempty(drx))
    drx=1;
end;
if (nargin<5)||(isempty(nswp))
    nswp=10;
end;
if (nargin<6)||(isempty(addswp))
    addswp=2;
end;

if (drx>m1)||(drx>m2)
    drx = min(m1,m2);
end;

rx=1;

curx = cell(2,1);
curx{1}=rand(m1,rx);
curx{2}=rand(m2,rx);
x = curx{1}*(curx{2}.');
derr = 2;
sp = 0;
resid_old = 1;
for swp=1:nswp
    

    % QR 2-1
    [q,rv]=qr(curx{2},0); % m2,rx' - rx',rx
    rx = size(q,2);
    curx{2}=q;
%     curx{1}=curx{1}*(rv.');
    
    % compute phi
    a2 = permute(mat{2}, [1 3 2]);
    a2 = reshape(a2, n2*ra, m2);
    phi = a2*curx{2}; % size n2*ra, rx
    phi = reshape(phi, n2, ra*rx);
    phi = (phi')*phi; % size ra*rx, ra*rx <-- for cplx should also work
    phi = reshape(phi, ra, rx, ra, rx);
    phi = reshape(permute(phi, [1 3 2 4]), ra*ra, rx*rx);
%     phi = reshape(permute(phi, [3 1 2 4]), ra, ra*rx*rx);
%     a2 = reshape(mat{2}, n2, m2*ra);
%     phi = (a2.')*phi; % size m2*ra, ra*rx
%     phi = reshape(phi, m2, ra*ra*rx);
%     phi = (curx{2}.')*phi; % size rx, ra*ra*rx
%     phi = reshape(permute(reshape(phi, rx, ra, ra, rx), [2 3 1 4]), ra*ra, rx*rx);
%     % And the projection of the matrix
    a1 = reshape(permute(mat{1}, [3 2 1]), ra*m1, n1);
    a1 = conj(a1)*reshape(mat{1}, n1, m1*ra); % size ra*m1, m1*ra <-- conjugate!
    a1 = reshape(a1, ra, m1, m1, ra);
    a1 = reshape(permute(a1, [1 4 2 3]), ra*ra, m1*m1);
%     a1 = reshape(permute(mat{1}, [3 2 1]), ra, m1*n1);
    a1 = (phi.')*a1; % size rx*rx, m1*n1
    a1 = reshape(a1, rx, rx, m1, m1);
    a1 = reshape(permute(a1, [2 4 1 3]), rx*m1, rx*m1);
%     a1 = reshape(mat{1}, n1*m1, ra)*phi; % size n1*m1, ra*rx*rx
%     a1 = reshape(a1, n1, m1, ra, rx, rx);
%     a1 = reshape(permute(a1, [1 3 2 4 5]), n1*ra, m1*rx*rx);
%     a1 = conj(reshape(permute(mat{1}, [2 1 3]), m1, n1*ra))*a1; % size m1, m1*rx*rx
%     a1 = reshape(a1, m1, m1, rx, rx);
%     a1 = reshape(permute(a1, [3 1 4 2]), rx*m1, rx*m1);
    
    %rhs: 
    
    rhs1 = rhs*conj(reshape(mat{2}, n2, m2*ra)); % size n1, m2*ra <-- conjugate
    rhs1 = reshape(rhs1, n1, m2, ra);
    rhs1 = reshape(permute(rhs1, [1 3 2]), n1*ra, m2);
    
    rhs1 = conj(reshape(permute(mat{1}, [2 1 3]), m1, n1*ra))*rhs1; % size m1, m2
%     rhs1 = rhs;
    rhs1 = rhs1*conj(curx{2}); % size m1,rx
    rhs1 = reshape(rhs1.', rx*m1, 1);
    
    curx{1}=a1 \ rhs1; % new first block
%     cond_a1 = cond(a1)
    curx{1}=reshape(curx{1}, rx, m1).';
    
    % Now, let's try the kickass by rank drx:
    if (mod(swp,addswp)==0)
%     if (sp>5)
        curx{1}=[curx{1}, randn(m1,drx)];
%         sp=0;
    end;
%     rx=rx+1;
    
    % Now, let's compute the second block
    [q,rv]=qr(curx{1},0); % m1,rx' - rx',rx
    rx = size(q,2);
    curx{1}=q;
    
    % compute phi
    a1 = permute(mat{1}, [1 3 2]);
    a1 = reshape(a1, n1*ra, m1);
    phi = a1*q; % size n1*ra, rx
    phi = reshape(phi, n1, ra*rx);
    phi = (phi')*phi; % size ra*rx, ra*rx
    phi = reshape(phi, ra, rx, ra, rx);
    phi = reshape(permute(phi, [1 3 2 4]), ra*ra, rx*rx);    
%     a1 = reshape(mat{1}, n1, m1*ra);
%     phi = (a1.')*phi; % size m1*ra, ra*rx
%     phi = reshape(phi, m1, ra*ra*rx);
%     phi = (curx{1}.')*phi; % size rx, ra*ra*rx
%     phi = reshape(permute(reshape(phi, rx, ra, ra, rx), [2 3 1 4]), ra*ra, rx*rx);
    % And the projection of the matrix
    a2 = reshape(permute(mat{2}, [3 2 1]), ra*m2, n2);
    a2 = conj(a2)*reshape(mat{2}, n2, m2*ra); % size ra*m2, m2*ra
    a2 = reshape(a2, ra, m2, m2, ra);
    a2 = reshape(permute(a2, [1 4 2 3]), ra*ra, m2*m2);
%     a2 = reshape(permute(mat{2}, [3 2 1]), ra, m2*n2);
    a2 = (phi.')*a2; % size rx*rx, m2*n2
    a2 = reshape(a2, rx, rx, m2, m2);
    a2 = reshape(permute(a2, [2 4 1 3]), rx*m2, rx*m2);
    
    %rhs: 
    rhs2 = rhs*conj(reshape(mat{2}, n2, m2*ra)); % size n1, m2*ra
    rhs2 = reshape(rhs2, n1, m2, ra);
    rhs2 = reshape(permute(rhs2, [1 3 2]), n1*ra, m2);
    
    rhs2 = conj(reshape(permute(mat{1}, [2 1 3]), m1, n1*ra))*rhs2; % size m1, m2
%     rhs2 = rhs;
    rhs2 = (curx{1}')*rhs2; % size rx,m2
    rhs2 = reshape(rhs2, rx*m2, 1);
    
    curx{2}=a2 \ rhs2; % new first block
    curx{2}=reshape(curx{2}, rx, m2).';
    
    x_new = curx{1}*(curx{2}.'); % size m1,m2
    derr = norm(x_new(:)-x(:))/norm(x(:));
    
    x = x_new;
    
    resid = norm(tt_mat_full_vec(mat, x(:))-rhs(:))/norm(rhs(:));
    conv_fact = resid_old/resid;
    if (conv_fact-1<1e-4)
        sp=sp+1;
    end;
    
    fprintf('als_solve: swp=%d, derr=%3.3e, rx=%d, resid=%3.3e, conv-1=%3.5e\n', swp, derr, rx, resid, conv_fact-1);
    if (derr<tol)
        break;
    end;
    resid_old = resid;
end;

x = x(:);

end
