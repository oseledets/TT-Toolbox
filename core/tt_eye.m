function [e]=tt_eye(n,d)
%[E]=TT_EYE(N,D)
% Identity matrix in TT format with mode size N^D
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


if (numel(n) == 1)
  n=n*ones(1,d);
end
e=cell(d,1);
for i=1:d
   e{i}=eye(n(i));
end
  return
end
