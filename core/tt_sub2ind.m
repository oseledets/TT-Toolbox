function [ind]=tt_sub2ind(siz, idx)
%Correct conversion of a multiindex to an index
%   [IND]=TT_SUB2IND(SIZ,IDX) 
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

d=numel(siz);
ind = 0;
shift=1;
for i=1:d
    ind=ind+shift*(idx(i)-1);
    shift = shift*siz(i);
end;
ind=ind+1;

end