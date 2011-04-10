function [newblock]=tt_addvecttobl(tt,block)
%[NEWBLOCK]=TT_ADDVECTTOBL(TT,BLOCK)
% Adds a tt-column to a tt-block
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

d=size(block,1);
if d==0
    d=size(tt,1);
    newblock=cell(d+1,1);
    for i=1:d-1
        newblock{i}=tt{i};
    end
    n=size(tt{d},1);
    r=size(tt{d},2);   
    newblock{d}=reshape(tt{d},n, r, 1);
    newblock{d+1}=ones(1,1);
else
    newblock=cell(d,1);
    
    newblock{1}=[block{1}, tt{1}];

    n=size(tt{d-1},1);
    rv2=size(tt{d-1},2);
    tt{d-1}=reshape(tt{d-1},[n,rv2,1]);

    for i=2:d-1
        n=size(tt{i},1);
        rv2=size(tt{i},2);
        rv3=size(tt{i},3);

        rb2=size(block{i},2);
        rb3=size(block{i},3);

        newblock{i}=zeros(n, rb2+rv2, rv3+rb3);
        newblock{i}(1:n,1:rb2,1:rb3)=block{i};
        newblock{i}(1:n,rb2+1:end,rb3+1:end)=tt{i};
    end
    
    nb=size(block{d},1);
    rb2=size(block{d},2);
    newblock{d}=[block{d}, zeros(nb,1); zeros(1,rb2), 1];

end

%    tt=tt_compr2(tt,eps);
return
end
