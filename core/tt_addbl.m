function [newblock]=tt_addbl(tt1,tt2)
%[NEWBLOCK]=TT_ADDBL(TT1,TT2)
% Concatenates two TT blocks
%
%
% TT Toolbox 1.1, 2009
%
%This is TT Toolbox, written by Ivan Oseledets, Olga Lebedeva
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------

d=size(tt1,1);
newblock=cell(d,1);
    
newblock{1}=[tt1{1}, tt2{1}];

 for i=2:d-1
        n=size(tt2{i},1);
        rv2=size(tt2{i},2);
        rv3=size(tt2{i},3);

        rb2=size(tt1{i},2);
        rb3=size(tt1{i},3);

        newblock{i}=zeros(n, rb2+rv2, rv3+rb3);
        newblock{i}(1:n,1:rb2,1:rb3)=tt1{i};
        newblock{i}(1:n,rb2+1:end,rb3+1:end)=tt2{i};
    end
    
    n1=size(tt1{d},1);
    r1=size(tt1{d},2); n2=size(tt2{d},1); r2=size(tt2{d},2);
    newblock{d}=[tt1{d}, zeros(n1,r2); zeros(n2,r1), tt2{d}];
return
end
