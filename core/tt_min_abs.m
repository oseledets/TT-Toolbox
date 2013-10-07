function [res,ind]=tt_min_abs(tt)
%Computes approximation to the minimal absolute value in a TT-tensor
%   [RES,IND]=TT_MIN_ABS(TT) Tries to compute the minimal absolute value in
%   TT-tensor
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

resmax = tt_max_abs(tt);
ons = tt_ones(size(tt));

[res,ind]=tt_max_abs(tt-ons*resmax);
res = resmax - res;

end
