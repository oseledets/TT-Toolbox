function [b]=chunk(b,i,j)
%Cut the (i,j) part out of the TT-tensor
%   [B]=CHUNK(B,I,J)
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