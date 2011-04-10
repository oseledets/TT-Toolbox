function [res] = tt_scal(tt,alpha)
%[RES]=TT_SCAL(TT,ALPHA)
%Multiply TT tensor by a scalar ALPHA
%
%
%
% TT Toolbox 1.0, 2009
%
%This is TT Toolbox, written by Ivan Oseledets.
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------
res=tt;
res{1}=res{1}*alpha;
end
