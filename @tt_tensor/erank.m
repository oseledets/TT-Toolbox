function [er]=erank(tt)
%Effective rank of the TT-tensor
%   [ER]=ERANK(TT)
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
%
% 31/01/2013 (Dmitry Savostyanov): 
% erank redefined like in tt90: erank(tt_ones)=1, etc.
%
n=tt.n; d=tt.d; r=tt.r;
if d<=1
    er=0.e0;
else
    sz=dot(n.*r(1:d),r(2:d+1));
    if sz == 0
        er=0.e0;
    else
        b=r(1)*n(1)+n(d)*r(d+1);
        if d == 2
            er=sz/b;
        else
            a=sum(n(2:d-1));
            er=(sqrt(b*b+4*a*sz)-b)/(2*a);
        end
    end
end
return
end