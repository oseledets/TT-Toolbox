function [res]=tt_mat_compr(mat,eps, max_rank)
%[RES]=TT_MAT_COMPR(MAT,EPS)
%Compress TT matrix with accuracy EPS
%
%
%
% TT Toolbox 1.1, 2009
%
%This is TT Toolbox, written by Ivan Oseledets.
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------
%n=size(mat{1},1);
%m=size(mat{1},2);

if (nargin<3)
    max_rank=[];
end;

d=size(mat,1);
n=zeros(d,1); m=zeros(d,1);
for i=1:d
  n(i)=size(mat{i},1); m(i)=size(mat{i},2);
end
res=tt_vec_to_mat(tt_compr2(tt_mat_to_vec(mat),eps,max_rank),n,m);
return
end
