function [sz] = size(tt)
%Mode sizes of the TT-matrix
%   [SZ]=SIZE(TT) Returns the mode sizes of the TT-matrix as a dx2 array of
%   integers
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
n=tt.n; m=tt.m;
n0=[n,m];
    sz=n0;
return
end