function y=tt_onlyone(idxs,arg)
%TT-Toolbox 2.2.2
% creates a tt_matrix with given core sizes arg 
%(arg is a 2-by-d array of core sizes or a different tt_matrix as an example)

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
