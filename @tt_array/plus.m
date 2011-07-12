function [a]=plus(b,c)
%[A]=PLUS(B,C)
%Adds two TT-array tensors B and C
if ( isempty(b) )
   a=c;
   return
elseif (isempty(c) )
   a=b;
   return
end

a=tt_array;
a.ns=b.ns;
% a.m=b.m;
% a.k=b.k;
a.tt=b.tt+c.tt;
return
end