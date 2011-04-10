function [res]=tt_kron(tt1,tt2)
%[RES]=TT_KRON(TT1,TT2)
% Kronecker product of two TT tensors, TT1,TT2
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
d1=size(tt1,1);
d2=size(tt2,1);
d=d1+d2;
res=cell(d,1);
for i=1:d1
    res{i}=tt1{i};
end
for i=1:d2
    res{i+d1}=tt2{i};
end
fx=res{d1+1};
ncur=size(fx,1);
r=size(fx,2);
res{d1+1}=reshape(fx,[ncur,1,r]);
return
end
