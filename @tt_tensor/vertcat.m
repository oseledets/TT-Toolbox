function [c]=vertcat(a,b)
%C=[A;B]
%   C=VERTCAT(A,B) Computes "addition of columns" of two TT-tensors A and B
%   The result has r(1) = sum(r(1,A)+r(1,B)). The last ranks should be
%   equal for A and B
%
%
% TT Toolbox 2.1, 2009-2012
%
%This is TT Toolbox, written by Ivan Oseledets et al.
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------


if (isempty(b) )
  c=a;
  return
elseif (isempty(a) )
  c=b;
  return
end
c=tt_tensor;
d=a.d;
ra=a.r;
n=a.n;
psa=a.ps;
cra=a.core;
rb=b.r;
crb=b.core;
psb=b.ps;
rc=ra+rb; 
if ( ra(d+1) == rb(d+1) )
rc(d+1)=ra(d+1); %The only modification
else
  error('Incorrect sizes in [A;B], please check');    
end
psc=cumsum([1;n.*rc(1:d).*rc(2:d+1)]);
crc=zeros(1,psc(d+1)-1);

%It seems: Ak
%          Bk --- the first core 
%all others are blocks

for i=1:d
  cr1=cra(psa(i):psa(i+1)-1);
  cr2=crb(psb(i):psb(i+1)-1);
  cr1=reshape(cr1,[ra(i),n(i),ra(i+1)]);
  cr2=reshape(cr2,[rb(i),n(i),rb(i+1)]);
  cr3=zeros(rc(i),n(i),rc(i+1));
  cr3(1:ra(i),:,1:ra(i+1))=cr1;
  cr3(rc(i)-rb(i)+1:rc(i),:,rc(i+1)-rb(i+1)+1:rc(i+1))=cr2;
  crc(psc(i):psc(i+1)-1)=cr3(:);
end
c.ps=psc;
c.d=d;
c.n=n;
c.r=rc;
c.core=crc;
