function [sz] = size(tt,p)
%Mode sizes of the TT-matrix
%   [SZ]=SIZE(TT) Returns the mode sizes of the TT-matrix as a dx2 array of
%   integers
%
%   [SZ]=SIZE(TT,P) Returns the mode-P sizes of the TT-matrix as an array
%   of d integers
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

if (nargin==1) || isempty(p)
    n=tt.n; m=tt.m;
    sz=[n,m];
else
    switch p
        case 1
            sz=tt.n;
        case 2
            sz=tt.m;
        otherwise
            error('tt_matrix/size : illegal value of parameter p.\n');
    end
end

return
end