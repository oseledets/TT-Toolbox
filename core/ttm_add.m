function [res]=ttm_add(mat1,mat2)
%[RES]=TTM_ADD(MAT1,MAT2)
%Adds two TTM matrices, MAT1, MAT2
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
n=tt_size(mat1);
m=tt_size(mat2);
res=tt_vec_to_mat(tt_add(tt_mat_to_vec(mat1),tt_mat_to_vec(mat2)),n,m);
return
end
