function [factors,cr]=tucker(tt,eps)
%Get Tucker factors and Tucker core in the TT-format
%   [FACTORS,CR]=TUCKER(TT,EPS) computes Tucker factors of a tensor
%   and the Tucker core in the TT-format with accuracy EPS
%
%
% TT Toolbox 2.1, 2009-2012
%
%This is TT Toolbox, written by Ivan Oseledets et al.
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------

%Zaglushka
[factors,cr]=tt_tuck(core(tt),eps);
cr=tt_tensor(cr);
return
end