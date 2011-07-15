function [c] = kron(a,b)
%[C]=KRON(A,B)
%Kronecker product of two TT-arrays
if ( isempty(a) )
  c=b;
  return;
elseif (isempty(b) )
  c=a;  
  return;
end
c=tt_array;
c.ns=[(a.ns);(b.ns)];
%c.d=a.d+b.d;
c.tt=kron(a.tt,b.tt);
return
end
