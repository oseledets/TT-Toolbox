function [x,y]=als_prec(mat,rhs,y,niter)
%Computes an approximate rank-1 solution to the linear system
%   [Y]=ALS_PREC(MAT,RHS,X) Computes (approximate) rank-1 solution to the
%   problems MAT*P = RHS, where A is a low-Kronecker rank matrix
%   P=X x Y, RHS is a low-Kronecker rank matrix
%   B is given as a two-2d cell array, the same is for X
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


a=mat{1}; %The first cell
b=mat{2}; 
R=size(a,3); %The rank
n=size(a,1);
m=size(b,1);
a1=permute(a,[1,3,2]); a1=reshape(a1,[n*R,n]);
b1=permute(b,[1,3,2]); b1=reshape(b1,[m*R,m]);
rhs1=rhs;
if ( nargin < 5 )
   niter=10;
end
u=rhs{1}; 
v=rhs{2}; 
rr=size(u,3);
u1=reshape(u,[n*n,rr]);
v1=reshape(v,[m*m,rr]);

for i=1:niter
    %fprintf('i=%d \n',i);
    %at=a1*x;
    %bt=b1*y;   
    %Iterate over x
    bt=b1*y;
    %bt is n*R*n, scalar product of 
    bt1=reshape(bt,[m,R,m]); bt1=permute(bt1,[1,3,2]); bt1=reshape(bt1,[m*m,R]);
    p=bt1'*bt1; %p matrix is ready
    %Compute the local matrix sum_{a,b} p_{ab} A^{\top}_b A_a,
    %it is A(k,i,b)*A(k,j,a)*p(a,b)
    a2=reshape(a,[n*n,R]); a2=a2*p; %a2 is (k,j,b), sum over (k,b) with (k,i,b)
    a3=reshape(a,[n,n,R]); a3=permute(a3,[2,1,3]); a3=reshape(a3,[n,n*R]);
    a2=reshape(a2,[n,n,R]); a2=permute(a2,[2,1,3]); a2=reshape(a2,[n,n*R]);
    loc_mat=a3*a2'; %Should be valid even for complex numbers
    %For the right hand side
    %Compute the q matrix <b,v> = b is m x m x R, V is m x m x rr
    %the result is R x rr
    q=bt1'*v1; %q is R x rr
    %Now compute right-hand side as
    %q(R,rr)*U(i,k,rr)*A(j,k,R)
    
    u2=u1*q'; % u2 is (i,k,R)
    u2=reshape(u2,[n,n*R]);
    a2=reshape(a,[n,n*R]);
    rhs=u2*a2';
    x= loc_mat \ rhs;
    %cond(loc_mat)
    %Iterate over y
    at=a1*x;
    %bt is n*R*n, scalar product of 
    at1=reshape(at,[n,R,n]); at1=permute(at1,[1,3,2]); at1=reshape(at1,[n*n,R]);
    p=at1'*at1; %p matrix is ready
    %Compute the local matrix sum_{a,b} p_{ab} A^{\top}_b A_a,
    %or the summation: \sum_{a,b,k} p(a,b) A(i,k,a)*A(j,k,b)
    b2=reshape(b,[m*m,R]); b2=b2*p; %a2 is (k,j,b), sum over (k,b) with (k,i,b)
    b3=reshape(b,[m,m,R]); b3=permute(b3,[2,1,3]); b3=reshape(b3,[m,m*R]);
    b2=reshape(b2,[m,m,R]); b2=permute(b2,[2,1,3]); b2=reshape(b2,[m,m*R]);
    loc_mat=b3*b2'; %Should be valid even for complex numbers
    %For the right hand side
    %Compute the q matrix <b,v> = b is m x m x R, V is m x m x rr
    %the result is R x rr
    q=at1'*u1; %q is R x rr
    %Now compute right-hand side as
    %q(R,rr)*U(i,k,rr)*A(j,k,R)
    
    v2=v1*q'; % u2 is (i,k,R)
    v2=reshape(v2,[m,m*R]);
    b2=reshape(b,[m,m*R]);
    rhs=v2*b2';
    y= loc_mat \ rhs;
    %    norm(tt_matrix(mat)*kron(tt_matrix(x,1e-10),tt_matrix(y,1e-10))-tt_matrix(rhs1))
    %keyboard;
    %norm(tt_matrix(mat)*kron(tt_matrix(x,1e-10),tt_matrix(y,1e-10))-tt_mat
    %rix(rhs1))
end
return
