function [res,ind]=tt_min(tt)
% Computes approximation to the minimal value in a TT-tensor
%   [RES,IND]=TT_MIN(TT) Tries to compute the minimal value in
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
[~, ind] = tt_max_abs(tt);
res = tt(ind);
bias = res;
while res > 0
	% Possibly we found maximum element instead of minimum.
	% Lets add a constant to the tensor and recalculate tt_max_abs.
	tt_biased = tt - bias;
	[~, ind] = tt_max_abs(tt_biased);
	res = tt_biased(ind);
	bias = bias + res;
end
res = tt(ind);
end