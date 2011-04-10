function [res]=tt_elem(tt,ind)
%[RES]=TT_ELEM(TT,IND)
%Returns the element of TT tensor in position given by multiindex IND
%
%
% TT Toolbox 1.1, 2009-2010
%
%This is TT Toolbox, written by Ivan Oseledets, Olga Lebedeva
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel
%
%For all questions, bugs and suggestions please mail
%ivan.oseledets@gmail.com
%---------------------------
d=size(tt,1);
mat=tt{1};
v0=mat(ind(1),:);
for k=2:d-1
 core=tt{k};
 core=squeeze(core); %To avoid rank-1 terms
 v0=v0*squeeze(core(ind(k),:,:));
end
mat=tt{d};
res=v0*mat(ind(d),:)';
return
end
