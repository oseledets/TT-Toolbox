function [tt]=row(tt,k)
%Returns the k-the row of a TT-tensor 
%   [TT]=ROW(TT,K) Returns the k-the row of a block-row TT-tensor
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


cr=tt.core;
ps=tt.ps;
r=tt.r;
n=tt.n;
cr1=cr(ps(1):ps(2)-1);
cr1=reshape(cr1,r(1),numel(cr1)/r(1));
cr1=cr1(k,:);
r(1)=size(cr1,1);
cr(ps(1):ps(2)-1)=[];
cr=[cr1(:); cr(:)];
%ps=[1;r(1)*n(1)*r(2);ps(3:end)-ps(2)+r(1)*n(1)*r(2)+1];
ps=cumsum([1;n.*r(1:end-1).*r(2:end)]);
tt.r=r;
tt.ps=ps;
tt.core=cr;
return
end
