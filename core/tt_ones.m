function [tt] = tt_ones(d,n)
%[TT]=TT_ONES(N,D)
% Tensor of all ones, modes size N, dimension D
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
 tt=cell(d,1);
 for k=1:d
    tt{k} = ones(n(k),1); %sqrt(n);
 end
return
end
