function [tt]=tt_transp(tt)
%Transposition of the TT matrix
%   [TT]=TT_TRANSP(TT) Transpose of the TT-matrix in TT1.0 format. Please 
%   avoid its usage: it will be removed in future releases. Use ' and .'
%   operators  from the object-oriented version
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
d=size(tt,1);
tt{1}=permute(tt{1},[2,1,3]);
tt{d}=permute(tt{d},[2,1,3]);
for k=2:d-1
  tt{k}=permute(tt{k},[2,1,3,4]);
end
return
end
