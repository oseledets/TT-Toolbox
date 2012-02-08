function [c] = kron(a,b)
%Kronecker product of two QTT_Tuckers
%   [C]=KRON(A,B) computes Kronecker product of A and B, where A and B are
%   in the QTT-Tucker format
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
