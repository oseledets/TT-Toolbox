function [mult]=sub_to_ind(ind,sz1)
%Converts a multiindex to a linear index
%   [MULT] = SUB_TO_IND(IND,SZ)
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
%---------------------------sz1=[1,cumprod(sz(1:numel(sz)-1))];
mult=dot((ind-1),sz1);%(ind-1).*sz;
mult = mult+1;
end