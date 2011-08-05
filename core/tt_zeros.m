function [tt] = tt_zeros(d,n)
%[TT]=TT_ZEROS(D,N)
%Creates a zero TT tensor of dimension D and mode size N
%
%
% TT Toolbox 1.1, 2009-2010
%
%This is TT Toolbox, written by Ivan Oseledets, Olga Lebedeva.
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------
if ( max(size(n)) == 1 )
   n=n*ones(1,d);
end
tt=cell(d,1);
tt{1}=zeros(n(1),1);
tt{d}=zeros(n(d),1);
for i=2:d-1
   tt{i}=zeros(n(i),1,1);
end
return
end
