function y = tt_blockwise(arr)
%% TT-Toolbox 2.2.1, 2009-2015
%This is TT Toolbox, written by Ivan Oseledets et al.
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel

%for questions on this particular code email Alexey Boyko, Skolkovo Institute of Science and Technology,
% a.boyko@skoltech.ru or alexey.boyko@skolkovotech.ru
%Moscow, Russia, 2016

%%
%INPUTS
% arr {2d cell array} of <tt_matrix>  variables
%OUTPUTS
% y <tt_matrix> 
%
%this function assembles a block matrix in tt-format from a given 2d cell array
%of blocks in tt-format

    [n,m]=size(arr); 

    block=zeros(n,m); block(1,1)=1;
    %since current version of TT-Toolbox doesn't support 
    %init like M=0; M=M+<tt_matrix>, one should construct 
    %a valid tt_matrix of consistent dimensions as an initial variable.
    
    
    y=tkron(arr{1,1},tt_matrix(block));

    for i=1:n
        for j=1:m
%             full(a{i,j})
            if (i~=1)||(j~=1)
                block=zeros(n,m); block(i,j)=1;
                y=y+tkron(arr{i,j},tt_matrix(block));
            end
        end
    end
return 
end