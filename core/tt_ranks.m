function [rks]=tt_ranks(tt)
%[RKS]=TT_RANKS(TT)
%Compute all ranks of TT decomposition
%Input : tensor TT
%Output: ranks RKS   
% 
%
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
d=size(tt,1);
if ( d == 1 ) 
  rks=[];
  return
end
rks=zeros(d-1,1);
rks(1)=size(tt{1},2);

for i=2:d-2
  rks(i) = size(tt{i},3);
end
rks(d-1)=size(tt{d},2);
return
end
