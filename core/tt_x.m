function [res]=tt_x(d,n)
%[RES]=TT_X(D,N)
% Function X in the QTT format, on a uniform grid with n^d points
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
%Function y=x in QTT format
res=cell(d,1);

res{1}=ones(n,2);
res{1}(:,1)=0:n-1;

shabl=zeros(n,2,2);
for j=1:n
    shabl(j,:,:)=eye(2,2);
end
ni=n;
for i=2:d-1
    res{i}=shabl;
    res{i}(:,2,1)=ni.*(0:n-1);
    ni=ni*n;
end

res{d}=ones(n,2);
res{d}(:,2)=ni.*(0:n-1);

return
end
