function [r]=rank(a,varargin)
%Computes ranks of the QTT-Tucker tensor
%   [R]=RANK(A) computes all ranks of the QTT-Tucker A
%
%   [R]=RANK(A,IND) computes the rank with index IND of the QTT-Tucker
%   tensor A
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

    r = 0;
    for i=1:a.dphys
        r(1:a.tuck{i}.d+1, i) = a.tuck{i}.r;
        if (i<a.dphys)
            r(a.tuck{i}.d+2, i) = a.core.r(i+1);
        end
    end
return
end