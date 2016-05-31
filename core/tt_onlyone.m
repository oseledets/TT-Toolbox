function y=tt_onlyone(idxs,arg)
% creates a tt_matrix with given core sizes arg 
%(arg is a 2-by-d array of core sizes or a different tt_matrix as an example)
%
%%TT-Toolbox 2.2.2, 2009-2016
%This is TT Toolbox, written by Ivan Oseledets et al.
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%for questions on this particular code email Alexey Boyko,
%Skolkovo Institute of Science and Technology,
%a.boyko@skoltech.ru or alexey.boyko@skolkovotech.ru
%Moscow, Russia, 2016

    if strcmp( class(arg),'tt_matrix')
        s=[arg.n arg.m].';
    else
        s=arg;
    end
    n=size(s,2);
    
    idx=indexify([idxs(1),idxs(2)],s);
    
    G={};
    for i=1:n
        Gbuf=zeros(1,s(1,i),s(2,i),1);
        Gbuf(1,idx(1,i),idx(2,i),1)=1;
        G{i}=Gbuf;
    end
    y=cell2core(tt_matrix,G);

end
