function [res]=full_to_nd(mat,eps,d) 
%[RES]=FULL_TO_ND(MAT,EPS,D)
%Converts 2^D x 2^D matrix into QTT format with accuracy EPS 
% TT Toolbox 1.0, 2009
%
%This is TT Toolbox, written by Ivan Oseledets.
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------
sz=2*ones(1,d*2);
res=reshape(mat,sz); %(i_1,...,i_d0,j_1,...,j_d0)
prm=zeros(1,d*2);
for i=1:d
  prm(2*i-1)=i;
  prm(2*i)=d+i;
end
res=permute(res,prm);
sz1=4*ones(1,d);
res=reshape(res,sz1);
res=full_to_tt(res,eps);
res=tt_vec_to_mat(res,2,2);
return
end
