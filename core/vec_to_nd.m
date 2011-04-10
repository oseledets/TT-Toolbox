function [tt] = vec_to_nd(vec,eps,d) 
%[TT]=VEC_TO_ND(VEC,EPS,D)
%Converts 2^d vector to TT format with accuracy EPS
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
tt=full_to_tt(reshape(vec,2*ones(1,d)),eps);
return
end
