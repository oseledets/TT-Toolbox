function [tt]=tt_transp(tt)
%[TT]=TT_TRANSP(TT)
%Transposition of the TT matrix
%
%
% TT Toolbox 1.1, 2009
%
%This is TT Toolbox, written by Ivan Oseledets, Olga Lebedeva
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
