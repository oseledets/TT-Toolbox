function [c] = power(a,b,varargin)
% C = A.^B
%   [C]=POWER(A,B) Compute A.^B for a TT-tensor A and natural B.
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

if isa(a,'tt_tensor') && numel(b) == 1
    c=a;
    b=b-1;
    while(b>0)
      c=c.*a;
      b=b-1;
    end
    
else
    error('Incorrect usage for TT-power');
end
