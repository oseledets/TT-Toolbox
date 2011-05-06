function [mat]=nd_to_full(tt)
%[MAT]=ND_TO_FULL(TT)
%Converts QTT representation (TT) to a 2^D x 2^D matrix
%Inverse of full_to_nd
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
ttm = tt_mat_to_vec(tt);
arr=tt_to_full(ttm);
sz=2*ones(1,d*2);
arr=reshape(arr,sz);
prm=zeros(1,d*2);
for i=1:d
  prm(2*i-1)=i;
  prm(2*i)=d+i;
end
mat=ipermute(arr,prm);
mat=reshape(mat,[2^d,2^d]);
return
end
