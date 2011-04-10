function [can_vec]=can_mat_to_vec(can_mat)
%[CAN_VEC]=CAN_MAT_TO_VEC(CAN_MAT)
%Converts canonical representation for matrix CAN_MAT to canonical
%representation of a long vector (CAN_VEC)
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
d=size(can_mat,1);
n=size(can_mat{1},1);
m=size(can_mat{1},2);
r=size(can_mat{1},3);
can_vec=cell(d,1);
for i=1:d
  can_vec{i} = reshape(can_mat{i},[n*m,r]);
end
return
end
