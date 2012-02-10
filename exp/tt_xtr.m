function [x]=tt_xtr(L,K,s)
%The X_s coordinate in the transposed QTT permutation of indices
%   [X]=TT_XTR(L,K,S)Computes QTT-representation of the vector X_s in the
%   transposed permutation with 2^L points in each mode 
%   and K variables. Ranks are still 2 here.
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

d=L*K;
r=2*ones(d+1,1);
r(1)=1;
r(d+1)=1;
n=2*ones(d,1);
ps=cumsum([1;n.*r(1:d).*r(2:d+1)]);
%cr=ones(ps(d+1)-1);
x0=zeros(1,2,2);
x0(:,1,:)=[0,1];
x0(:,2,:)=[0,1];
cr(ps(1):ps(2)-1)=x0(:);
x0=zeros(2,2,1);
x0(:,1,:)=[1;0];
x0(:,2,:)=[1;0];
cr(ps(d):ps(d+1)-1)=x0(:);
for i=2:d-1
   x0=zeros(2,2,2);
   x0(:,1,:)=eye(2);
   x0(:,2,:)=eye(2);
   cr(ps(i):ps(i+1)-1)=x0(:);
end
x=tt_tensor;

%Now distribute them in right places: the cores in the 
%transposed permutation are with numbers K*(i-1)+s, i=1,...,L

for i=1:L
   curI=K*(i-1)+s;
   if ( curI ~= 1 && curI ~= d )
      x0=zeros(2,2,2);
      x0(:,1,:)=[1,0;0,1];
      x0(:,2,:)=[1,0;2^(i-1),1];
   elseif ( curI == 1 )
      x0=zeros(1,2,2);
      x0(:,1,:)=[0,1];
      x0(:,2,:)=[2^(i-1),1];
   else
      x0=zeros(2,2,1);
      x0(:,1,:)=[1;0];
      x0(:,2,:)=[1;2^(i-1)]; 
   end
   cr(ps(curI):ps(curI+1)-1)=x0(:);
end



x.d=d;
x.n=n;
x.ps=ps;
x.core=cr;
x.r=r;
end