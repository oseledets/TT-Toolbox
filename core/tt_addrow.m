function [newblock]=tt_addrow(tt,block)
%[NEWBLOCK]=TT_ADDROW(TT,BLOCK)
% Adds a tt-row to a tt-blockrow
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

d=size(block,1);
if d==0
    d=size(tt,1);
    newblock=cell(d+1,1);
    for i=2:d+1
        newblock{i}=tt{i-1};
    end
    n=size(tt{1},1);
    r=size(tt{1},2);   
    newblock{2}=reshape(tt{1},[n,1,r]);
    newblock{1}=ones(1,1);
else
    newblock=cell(d,1);
    
    newblock{d}=[block{d}, tt{d-1}];

    n=size(tt{1},1);
    rv2=size(tt{1},2);
    tt{1}=reshape(tt{1},[n,1,rv2]);

    for i=2:d-1
        n=size(tt{i-1},1);
        rv2=size(tt{i-1},2);
        rv3=size(tt{i-1},3);

        rb2=size(block{i},2);
        rb3=size(block{i},3);

        newblock{i}=zeros(n, rb2+rv2, rv3+rb3);
        newblock{i}(1:n,1:rb2,1:rb3)=block{i};
        newblock{i}(1:n,rb2+1:end,rb3+1:end)=tt{i-1};
    end
    
    nb=size(block{1},1);
    rb2=size(block{1},2);
    newblock{1}=[block{1}, zeros(nb,1); zeros(1,rb2), 1];

end

return
end
