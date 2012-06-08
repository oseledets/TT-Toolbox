function [c] = kron(a,b,typ)
%Kronecker product of two TT-matrices
%   [C]=KRON(A,B,TYP) Kronecker product of two TT-matrices. One of the
%   arguments can be empty. In this case, the other, nonempty argument is
%   returned. TYP can be absent, then it is old-tt-style Kronecker product 
%   (not the standard definition). If TYP is present, then it is 
%   consistent with the standard definition. This may change 
%   change in future to use only the standard definition
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
  return;
elseif (isempty(b) )
  c=a;  
  return;
end
if ( nargin < 3 )
    c=tt_matrix;
    c.n=[(a.n);(b.n)];
    c.m=[(a.m);(b.m)];
    c.tt=kron(a.tt,b.tt);
else
   c = kron(b,a);
end
return
end
