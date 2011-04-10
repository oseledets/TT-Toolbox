function [tt]=tt_lapl(d)
%[TT]=TT_LAPL(D)
%Compute a standard second order difference
%of size 2^D x 2^D in TT format
%
%
% TT Toolbox 1.1, 2009-2010
%
%This is TT Toolbox, written by Ivan Oseledets, Olga Lebedeva.
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------
tt1=tt_shf(d); tt2=tt_transp(tt1);
tt1=tt_mat_to_vec(tt1); tt2=tt_mat_to_vec(tt2);
tt2=tt_add(tt1,tt2);
tt=tt_add(tt_scal(tt_mat_to_vec(tt_eye(2,d)),2),tt_scal(tt2,-1.0));
tt=tt_compr2(tt,1e-13);
tt=tt_vec_to_mat(tt,2,2);
return
end
