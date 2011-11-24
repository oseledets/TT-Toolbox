function [c]=horzcat(a,b)
%C=HORZCAT(A,B)
%C=[A,B]
%Computes "addition of rows" of two TT-tensors A and B
%The result has r(d+1) = sum(r(d+1,A)+r(d+1,B))
%r(1) should be equal to the same of A & B

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
if ( ra(1) == rb(1) )
rc(1)=ra(1); %The only modification
else
  error('Incorrect sizes in [A,B], please check');    
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
