function [mat]=nd_to_full(tt)
%[MAT]=ND_TO_FULL(TT)
%Converts QTT representation (TT) to a 2^D1 x 2^D2 matrix
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
n1 = ones(1,d); n2 = ones(1,d); n = ones(1,2*d);
for i=1:d
    n1(i) = size(tt{i},1);    
    n2(i) = size(tt{i},2);
    n(2*i-1)=n1(i); n(2*i)=n2(i);
end;
ttm = tt_mat_to_vec(tt);
arr=tt_to_full(ttm);
% sz=2*ones(1,d*2);
arr=reshape(arr,n);
prm=zeros(1,d*2);
for i=1:d
  prm(2*i-1)=i;
  prm(2*i)=d+i;
end
mat=ipermute(arr,prm);
mat=reshape(mat,[prod(n1),prod(n2)]);
return
end
