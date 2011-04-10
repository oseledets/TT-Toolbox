function [tt]=tt_random(n,d,r)
%[TT]=TT_RANDOM(N,D,R)
%Generates a random D-dimensional tensor TT with mode sizes N and ranks R
%N can be either mode vector or a number
%R can be either mode vector or a number
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
if ( max(size(n)) == 1 )
  n = n*ones(1,d);
end
if ( max(size(r)) == 1 )
  r = r*ones(1,d-1);
end
tt=cell(d,1);
tt{1}=randn(n(1),r(1));
tt{d}=randn(n(d),r(d-1));
for k=2:d-1
    tt{k}=randn(n(k),r(k-1),r(k));
end
return
end
