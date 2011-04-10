function [tt]=tt_column(block,i)
%Takes a column from the tt-block
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
d=size(block,1)-1;

tt=cell(d,1);
for k=1:d-1
    tt{k}=block{k};
end
rv1=size(block{d},1);
rv2=size(block{d},2);
rv3=size(block{d},3);
tt{d}=reshape(block{d},[rv1*rv2,rv3])*block{d+1}(i,:)';
tt{d}=reshape(tt{d},[rv1,rv2]);

%    tt=tt_compr2(tt,eps);
return
end
