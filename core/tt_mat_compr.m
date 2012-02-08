function [res]=tt_mat_compr(mat,eps, max_rank)
%Tensor rounding for the TT-matrix in TT1.0 format
%   [RES]=TT_MAT_COMPR(MAT,EPS) Compress TT matrix with accuracy EPS. Avoid
%   the usage, use round() function of the object-oriented design. Will be
%   removed in future releases.
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
