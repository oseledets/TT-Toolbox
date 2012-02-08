

function [M]=Fd_mtx(d, a, bound)
%Finite difference approximation of scalar diffusion equation in QTT
%   [M]=FD_MTX(D, A, BOUND)
%   Generate finite difference matrix from diffusion coefficient A(n,n) or A(n,n,n)
%   M - n^D-by-n^D sparse matrix
%   D - dimensions, 
%   BOUND:
%       0 - Dirichlet,
%       1 - Neuman

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


n = size(a, 1)-2;
if (d==1) 
  bar_a = zeros(n,3);
  for i=1:n
     
     bar_a(i,2)=(a(i+2)+2*a(i+1)+a(i))*0.25;
     bar_a(i,1)=(a(i)+a(i+1))*0.5;
     bar_a(i,3)=(a(i+1)+a(i+2))*0.5;
  end
  D1=2*eye(n)-diag(ones(n-1,1),1)-diag(ones(n-1,1),-1);
  %M=D1*(n+1)^2;
  M=D1;
  bar_a_M=spdiags(bar_a,[-1,0,1],n,n);
  M = M.*bar_a_M;
  return;
end

if (d==2)
    if (bound==0)
        bar_a = zeros(n^2, 5);
        for i=1:n
            for j=1:n
                bar_a((i-1)*n+j, 3) = (a(i+1,j+1)+a(i+2,j+1)+a(i+1,j+2)+a(i+2,j+2))*0.25;
                if (j>1) bar_a((i-1)*n+j-1, 2) = (a(i+1,j+1)+a(i+2,j+1))*0.5; end;
                if (j<n) bar_a((i-1)*n+j+1, 4) = (a(i+1,j+2)+a(i+2,j+2))*0.5; end;
                if (i>1) bar_a((i-1-1)*n+j, 1) = (a(i+1,j+1)+a(i+1,j+2))*0.5; end;
                if (i<n) bar_a((i+1-1)*n+j, 5) = (a(i+2,j+1)+a(i+2,j+2))*0.5; end;
            end;
        end;
        D1 = spdiags(-1*ones(n,1), [-1], n, n)+spdiags(2*ones(n,1), [0], n, n)+spdiags(-1*ones(n,1), [1], n, n);
        M = kron(D1*(n+1)^2, speye(n))+kron(speye(n), D1*(n+1)^2);
        bar_a_M = spdiags(bar_a, [-n, -1, 0, 1, n], n^2, n^2);
        M = M.*bar_a_M;
    else
        wa = zeros(n+3,n+3);
        wa(1:n+2,1:n+2)=a;
        wa(n+3,1:n+2)=a(n+2,:);
        wa(1:n+2,n+3)=a(:,n+2);
        wa(n+3,n+3)=a(n+2,n+2);
        n = n+2;
        bar_a = zeros(n^2, 5);
        for i=1:n
            for j=1:n
                bar_a((i-1)*n+j, 3) = (wa(i,j)+wa(i+1,j)+wa(i,j+1)+wa(i+1,j+1))*0.25;
                if (j>1) bar_a((i-1)*n+j-1, 2) = (wa(i,j)+wa(i+1,j))*0.5; end;
                if (j<n) bar_a((i-1)*n+j+1, 4) = (wa(i,j+1)+wa(i+1,j+1))*0.5; end;
                if (i>1) bar_a((i-1-1)*n+j, 1) = (wa(i,j)+wa(i,j+1))*0.5; end;
                if (i<n) bar_a((i+1-1)*n+j, 5) = (wa(i+1,j)+wa(i+1,j+1))*0.5; end;
            end;
        end;        
        D1 = spdiags(-1*ones(n,1), [-1], n, n)+spdiags(2*ones(n,1), [0], n, n)+spdiags(-1*ones(n,1), [1], n, n);
        D1(1,1)=1;
        D1(n,n)=1;
        M = kron(D1*(n-1)^2, speye(n))+kron(speye(n), D1*(n-1)^2);
        bar_a_M = spdiags(bar_a, [-n, -1, 0, 1, n], n^2, n^2);
        M = M.*bar_a_M;
    end;
end;

if (d==3)
    if (bound==0)
        bar_a = zeros(n^3, 7);
        wa=a;
        for i=1:n
            for j=1:n
                for k=1:n
                    bar_a((i-1)*n^2+(j-1)*n+k, 4) = (wa(i+1,j+1,k+1)+wa(i+2,j+1,k+1)+wa(i+1,j+2,k+1)+wa(i+2,j+2,k+1)+wa(i+1,j+1,k+2)+wa(i+2,j+1,k+2)+wa(i+1,j+2,k+2)+wa(i+2,j+2,k+2))*0.125;
                    if (k>1) bar_a((i-1)*n^2+(j-1)*n+k-1, 3) = (wa(i+1,j+1,k+1)+wa(i+2,j+1,k+1)+wa(i+1,j+2,k+1)+wa(i+2,j+2,k+1))*0.25; end;
                    if (k<n) bar_a((i-1)*n^2+(j-1)*n+k+1, 5) = (wa(i+1,j+1,k+2)+wa(i+2,j+1,k+2)+wa(i+1,j+2,k+2)+wa(i+2,j+2,k+2))*0.25; end;
                    
                    if (j>1) bar_a((i-1)*n^2+(j-1-1)*n+k, 2) = (wa(i+1,j+1,k+1)+wa(i+2,j+1,k+1)+wa(i+1,j+1,k+2)+wa(i+2,j+1,k+2))*0.25; end;
                    if (j<n) bar_a((i-1)*n^2+(j+1-1)*n+k, 6) = (wa(i+1,j+2,k+1)+wa(i+2,j+2,k+1)+wa(i+1,j+2,k+2)+wa(i+2,j+2,k+2))*0.25; end;
                    
                    if (i>1) bar_a((i-1-1)*n^2+(j-1)*n+k, 1) = (wa(i+1,j+1,k+1)+wa(i+1,j+2,k+1)+wa(i+1,j+1,k+2)+wa(i+1,j+2,k+2))*0.25; end;
                    if (i<n) bar_a((i+1-1)*n^2+(j-1)*n+k, 7) = (wa(i+2,j+1,k+1)+wa(i+2,j+2,k+1)+wa(i+2,j+1,k+2)+wa(i+2,j+2,k+2))*0.25; end;                    
                end;
            end;
        end;        
        D1 = spdiags(-1*ones(n,1), [-1], n, n)+spdiags(2*ones(n,1), [0], n, n)+spdiags(-1*ones(n,1), [1], n, n);
        M = kron(kron(D1*(n+1)^2, speye(n)), speye(n))+kron(kron(speye(n), D1*(n+1)^2), speye(n)) + kron(kron(speye(n), speye(n)), D1*(n+1)^2);
        bar_a_M = spdiags(bar_a, [-n^2, -n, -1, 0, 1, n, n^2], n^3, n^3);
        M = M.*bar_a_M;
    else
        wa = zeros(n+3,n+3,n+3);
        wa(1:n+2, 1:n+2, 1:n+2) = a;
        wa(n+3, 1:n+2, 1:n+2) = a(n+2, :, :);
        wa(1:n+2, n+3, 1:n+2) = a(:, n+2, :);
        wa(1:n+2, 1:n+2, n+3) = a(:, :, n+2);
        wa(1:n+2, n+3, n+3) = a(:, n+2, n+2);
        wa(n+3, 1:n+2, n+3) = a(n+2, :, n+2);
        wa(n+3, n+3, 1:n+2) = a(n+2, n+2, :);
        wa(n+3,n+3,n+3)=a(n+2,n+2,n+2);
        n = n+2;
        bar_a = zeros(n^3, 7);
        for i=1:n
            for j=1:n
                for k=1:n
                    bar_a((i-1)*n^2+(j-1)*n+k, 4) = (wa(i,j,k)+wa(i+1,j,k)+wa(i,j+1,k)+wa(i+1,j+1,k)+wa(i,j,k+1)+wa(i+1,j,k+1)+wa(i,j+1,k+1)+wa(i+1,j+1,k+1))*0.125;
                    if (k>1) bar_a((i-1)*n^2+(j-1)*n+k-1, 3) = (wa(i,j,k)+wa(i+1,j,k)+wa(i,j+1,k)+wa(i+1,j+1,k))*0.25; end;
                    if (k<n) bar_a((i-1)*n^2+(j-1)*n+k+1, 5) = (wa(i,j,k+1)+wa(i+1,j,k+1)+wa(i,j+1,k+1)+wa(i+1,j+1,k+1))*0.25; end;
                    
                    if (j>1) bar_a((i-1)*n^2+(j-1-1)*n+k, 2) = (wa(i,j,k)+wa(i+1,j,k)+wa(i,j,k+1)+wa(i+1,j,k+1))*0.25; end;
                    if (j<n) bar_a((i-1)*n^2+(j+1-1)*n+k, 6) = (wa(i,j+1,k)+wa(i+1,j+1,k)+wa(i,j+1,k+1)+wa(i+1,j+1,k+1))*0.25; end;
                    
                    if (i>1) bar_a((i-1-1)*n^2+(j-1)*n+k, 1) = (wa(i,j,k)+wa(i,j+1,k)+wa(i,j,k+1)+wa(i,j+1,k+1))*0.25; end;
                    if (i<n) bar_a((i+1-1)*n^2+(j-1)*n+k, 7) = (wa(i+1,j,k)+wa(i+1,j+1,k)+wa(i+1,j,k+1)+wa(i+1,j+1,k+1))*0.25; end;                    
                end;
            end;
        end;        
        D1 = spdiags(-1*ones(n,1), [-1], n, n)+spdiags(2*ones(n,1), [0], n, n)+spdiags(-1*ones(n,1), [1], n, n);
        D1(1,1)=1;
        D1(n,n)=1;
        M = kron(kron(D1*(n-1)^2, speye(n)), speye(n))+kron(kron(speye(n), D1*(n-1)^2), speye(n)) + kron(kron(speye(n), speye(n)), D1*(n-1)^2);
        bar_a_M = spdiags(bar_a, [-n^2, -n, -1, 0, 1, n, n^2], n^3, n^3);
        M = M.*bar_a_M;
    end;    
end;

end