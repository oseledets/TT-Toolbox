function [c] = kron(a,b)
%[C]=KRON(A,B)
%Kronecker product of two QTT_tuckers

%Time to replace it by something meaningful
if ( isempty(a) )
  c=b;
  return
elseif ( isempty(b) )
  c=a;
  return
end

c = qtt_tucker;
c.dphys = a.dphys+b.dphys;
c.tuck = cell(c.dphys, 1);
c.core = kron(a.core, b.core);
c.tuck(1:a.dphys) = a.tuck;
c.tuck(a.dphys+1:c.dphys) = b.tuck;

end
