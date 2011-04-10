function [tt] = can_to_tt(can)
%[TT]=CAN_TO_TT(CAN)
%Converts the canonical representation of a tensor CAN
%to TT representation
%No compression is done; use TT_COMPR2 to compress
%the result
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
d=size(can,1);
tt=cell(d,1);
tt{1}=can{1};
r=size(can{1},2);
tt{d}=can{d}; % Nothing more can be done 
for i=2:d-1
    fact=can{i};
    ncur=size(fact,1); %Matrix is ncur x r
    core=zeros(ncur,r,r);
    for s=1:r
       for j=1:ncur
           core(j,s,s)=fact(j,s);
       end 
    end 
    tt{i}=core;
end
return
end
