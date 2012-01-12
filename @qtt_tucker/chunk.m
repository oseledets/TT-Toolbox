function [b]=chunk(b,i,j)
%[B]=CHUNK(B,I,J)
%Cut the (i,j) part out of the TT-tensor
if ( i > j ) 
 tmp=j; 
 j=i;
 i=tmp;
end
ps=b.ps;
r=b.r;
n=b.n;
cr=b.core;
r=r(i:j+1); 
n=n(i:j);
cr=cr(ps(i):ps(j+1)-1);
ps=ps(i:j+1);
ps=ps-ps(1)+1;
b.r=r;
b.n=n;
b.ps=ps;
b.core=cr;
b.d=numel(n);
return
end