function [res]=can_eye(n,d)
%[RES]=CAN_EYE(N,D)
%Identity matrix in the canonical format
%
%
% TT Toolbox 1.1, 2009-2010
%
%This is TT Toolbox, written by Ivan Oseledets, Olga Lebedeva
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------
if ( numel(n) == 1 )
   n=n*ones(1,d);   
end
res=cell(d,1);

for i=1:d
    res{i}=eye(n(i));
end
return
end
