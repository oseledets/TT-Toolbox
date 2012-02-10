function [res] = tt_scal(tt,alpha)
%Multiply tensor by a scalar in TT1.0 format
%   [RES]=TT_SCAL(TT,ALPHA) Multiply TT tensor by a scalar ALPHA in TT1.0
%   format. Please avoid its usage: it will be removed in
%   future releases. Use * operator from the object-oriented version
%
%
% TT-Toolbox 2.2, 2009-2012
%
%This is TT Toolbox, written by Ivan Oseledets et al.
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please   
%ivan.oseledets@gmail.com
%---------------------------
res=tt;
res{1}=res{1}*alpha;
end
