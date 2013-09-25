function [tt]=kron2(tt1,tt2)
%Kronecker product of two TT-matrices in non-standard ordering
%   [TT]=KRON2(TT1,TT2) Computes Kronecker product of two TT-matrices in
%   non-standard ordering: row indices of TT1 and TT2 are "mixed". The
%   ranks are multiplied after this operation. It can be useful for the
%   matrix inversion
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

if ( isempty(tt1) )
  tt=tt2;
  return
elseif ( isempty(tt2) )
  tt=tt1;
  return
end

n1=tt1.n; m1=tt1.m; n2=tt2.n; m2=tt2.m;
t1=tt1.tt; t2=tt2.tt;
cr1=t1.core; cr2=t2.core; r1=t1.r; r2=t2.r; ps1=t1.ps; ps2=t2.ps;
d=t1.d;
r=r1.*r2;
n=n1.*n2;  
m=m1.*m2;
sz=dot(n.*r(1:d),r(2:d+1));
p=n.*m;
  pos=cumsum([1;p.*r(1:d).*r(2:d+1)]);

cr=zeros(sz,1);
for i=1:d
   core1=cr1(ps1(i):ps1(i+1)-1);
   core2=cr2(ps2(i):ps2(i+1)-1);
   core1=core1(:);
   core2=core2(:);
   cr0=core1*core2.';
   cr0=reshape(cr0,[r1(i),n1(i),m1(i),r1(i+1),r2(i),n2(i),m2(i),r2(i+1)]);
   cr0=permute(cr0,[1,5,2,6,3,7,4,8]);
   cr(pos(i):pos(i+1)-1)=cr0(:);
end
t=tt_tensor;
t.core=cr;
t.ps=pos;
t.r=r;
t.n=p;
t.d=d;
tt=tt_matrix;
tt.tt=t;
tt.n=n;
tt.m=m;
return
end
