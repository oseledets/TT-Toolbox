function x_ = permute(x_, order, eps)
%Permute dimensions.
%   Y = permute(X, ORDER, EPS) permutes the dimensions of the TT-tensor X
%   according to ORDER, delivering a result at relative accuracy EPS. This
%   function is equivalent to 
%      Y = tt_tensor(permute(reshape(full(X), X.n'),ORDER), EPS) 
%   but avoids the conversion to the full format.
%
%
% Simon Etter, Summer 2015
% Seminar of Applied Mathematics, ETH Zurich
%
%
% TT-Toolbox 2.2, 2009-2012
%
%This is TT Toolbox, written by Ivan Oseledets et al. Institute of
%Numerical Mathematics, Moscow, Russia webpage:
%http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com 
%---------------------------
x_.tt = permute(x_.tt, order, eps);
x_.n = x_.n(order);
x_.m = x_.m(order);
end
