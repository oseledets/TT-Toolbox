function [rks]=tt_ranks(tt)
%Compute all ranks of the TT-decomposition in TT1.0 format
%   [RKS]=TT_RANKS(TT) Computes all ranks of TT decomposition in TT1.0
%   format. Please avoid its usage: it will be removed in
%   future releases. Use rank() from the object-oriented version.
%
%
% TT-Toolbox 2.2, 2009-2012
%
%This is TT Toolbox, written by Ivan Oseledets et al.
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
