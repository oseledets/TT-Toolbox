function [can_res] = can_add(can1,can2)
%
% TT Toolbox 1.0, 2009
%
%This is TT Toolbox, written by Ivan Oseledets.
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------

d=size(can1,1);
can_res=cell(d,1);
for i=1:d
    can_res{i}=[can1{i} can2{i}];
end
return
end
