function [a]=plus(b,c)
%[A]=PLUS(B,C)
%Adds two QTT_TUCKER Tensors B and C
if ( isempty(b) )
 a=c;
 return
elseif ( isempty(c) )
 a=b;
 return
end
a=qtt_tucker;
cb=b.core; cc=c.core;
rb=rank(cb); rc=rank(cc);
nb=size(cb); nc=size(cc);
tb=b.tuck; tc=c.tuck;
d=b.dphys;
%Determine size of the result
ra=rb+rc;
na=nb+nc;
if ( rc(1) == rb(1) )
    ra(1)=rb(1);
else
  error('Inconsistent sizes in mode 1');    
end
if ( rc(d+1) == rb(d+1) )
  ra(d+1)=rb(d+1);
else
  error('Inconsistent sizes in the last mode');    
end
corea=[];
for i=1:d
  c1=cb{i}; c2=cc{i};
  c1=reshape(c1,[rb(i),nb(i),rb(i+1)]);
  c2=reshape(c2,[rc(i),nc(i),rc(i+1)]);
  cres=zeros(ra(i),na(i),ra(i+1));
  cres(1:rb(i),1:nb(i),1:rb(i+1))=c1;
  cres((ra(i)-rc(i)+1):ra(i),(na(i)-nc(i)+1):na(i),(ra(i+1)-rc(i+1)+1):ra(i+1))=c2;
  corea=[corea; cres(:)];
end
ca=tt_tensor; 
ca.d=d;
ca.n=na;
ca.core=corea;
ca.r=ra;
ca.ps=cumsum([1;ca.n.*ca.r(1:ca.d).*ca.r(2:ca.d+1)]);
for i=1:d
   ta{i}=[tb{i},tc{i}];
end
a.core=ca;
a.tuck=ta;
a.dphys=d;
a.sz=b.sz;
return
end