function [mem]=tt_mem(tt)
%[MEM]=TT_MEM(TT)
% Memory required to store the tensor TT
%
%
% TT Toolbox 1.0, 2009-2010
%
%This is TT Toolbox, written by Ivan Oseledets, Olga Lebedeva
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------
mem=0;
d=size(tt,1);
for k=1:d
    mem=mem+numel(tt{k});
end
end
