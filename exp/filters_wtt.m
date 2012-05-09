function [wtt_transform] = filters_wtt (a,rmax,eps,sz)
%Computes the WTT transform filters for a given tensor
%   [WTT_TRANSFORM]=FILTERS_WTT(A,RMAX,EPS,SZ)
%   Computes the set of filters for transform of A, A should be given 
%   either as d-dimensional array, or as by size-vector SZ
%   WTT_TRANSFORM is a structure with two fields: 
%   u and rks
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

if ( nargin == 1 )
sz=size(a); d = numel(sz);
else
   d = numel(sz);
end
N=numel(a);
rks=zeros(d,1); u=cell(d,1);
rks(1)=1;
a=reshape(a,[sz(1),N/sz(1)]);
for k=1:d
   [m,p]=size(a);
   if ( m > p ) 
      [u1,s1,v1]=svd(a);
   else
      [u1,s1,v1]=svd(a,'econ');
   end
   
   %Fix for u1 SLIGHTLY
%   a1=a;
%    if ( m > p )
%  for s=1:5
%    b1=u1'*a1;
%    m1=size(b1,2);
%    sp=zeros(m,1);
%    for q=1:m1
%      b=b1(:,q); 
%      sp(q)=norm(b,2)/norm(b,1); %this is a sparsity measure: for one 1 it is 1, for all ones --- 1/sqrt(np)
%    end
%    sgm=0.95; %cutoff for sparse vectors
%    ind1=find(sp>sgm);
%    if ( numel(ind1) > 0 )
%     a1=a1(:,ind1);
%     if ( m > p ) 
%       [u1,s1,v1]=svd(a1);
%     else
%       [u1,s1,v1]=svd(a1,'econ');
%     end
%    end
%  end
%    end
    s1=diag(s1);
   %s1
   r0=my_chop2(s1,eps*norm(s1)/sqrt(d-1));
   r=min([m,p,rmax,r0]);
   N = N*r/(sz(k)*rks(k));
   rks(k+1)=r;
   a=u1'*a; 
   if ( k < d )
    a=a(1:r,:); a=reshape(a,[r*sz(k+1),N/(r*sz(k+1))]);
   end
   u{k}=u1;
   % fprintf('k=%d r=%d Size of u = %d x %d, Size of a is %d x %d \n',k,r,size(u1,1),size(u1,2),size(a,1),size(a,2));
end
rks(d+1)=1;
%u{d}=v1;
wtt_transform=struct;
wtt_transform.u=u;
wtt_transform.rks=rks;
wtt_transform.sz=sz;
return
end
