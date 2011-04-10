function [res]=tt_dist2(tt1,tt2)
%[RES]=TT_DIST2(TT1,TT2)
%More accurate and but more expensive than
%TT_DIST computation of a distance between two tensors
% TT1 and TT2 in Frobenius norm: just recompression
% of their difference
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
rs=tt_add(tt1,tt_scal(tt2,-1.0));
rs=tt_compr2(rs,0.d0); %After recompression cores are 
d=size(rs,1);
res=norm(rs{1}(:))/sqrt(size(rs{1},2));
res=res*norm(rs{d}(:));
for k=2:d-1
  res = res * norm(rs{k}(:))/sqrt(size(rs{k},3));
end
return
end
