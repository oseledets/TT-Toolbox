function [res,ind]=tt_max(tt)
%Compute the (supposedly) maximal element in a TT-tensor
%   [RES,IND]=TT_MAX(TT) Compute element quasi-maximal in 
%   a TT tensor. It return the value RES and the multiindex
%   IND where this element is located.
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
[~, ind] = tt_max_abs(tt);
res = tt(ind);
if res < 0
	% Possibly we found minimum element instead of maximum.
	% Lets add a constant to the tensor and recalculate tt_max_abs.
	tt_tmp = tt - res;
	[~, ind] = tt_max_abs(tt_tmp);
	res = tt(ind);
end
end