function [mult]=sub_to_ind(ind,sz)
%[MULT] = SUB_TO_IND(IND,SZ)
%Converts a multiindex to a linear index
% TT Toolbox 1.1, 2009-2011
%
%This is TT Toolbox, written by Ivan Oseledets, Olga Lebedeva
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------
sz1=[1,cumprod(sz(1:numel(sz)-1))];
mult=dot((ind-1),sz1);%(ind-1).*sz;
mult = mult+1;
%d=size(sz,2);
%prs=1;
%mult=1;
%for k=1:d
%    mult=mult+(ind(k)-1)*prs;
%    prs = prs*sz(k);
%end
return
