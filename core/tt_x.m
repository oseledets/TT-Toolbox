function [res]=tt_x(n,varargin)
%Computes X in the QTT-format
%   [RES]=TT_X(N,D), or [RES]=TT_X(N): Vector with elements 
%   0:N(1)*...*N(D)-1 is transformed into a tensor of size 
%   N(1) x N(2) x ... x N(D)
%   This function returns its TT-decomposition with rank 2
%   This function is useful for the computation of the QTT-decomposition
%   of the function y=x
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
if (numel(n) == 1)
    if (numel(varargin)>0)
        d=varargin{1};
    else
        d=1;
    end;
    n=n*ones(1,d);
else
  d=numel(n);
end



res=cell(d,1);

ni=n(1);

for i=2:d-1
    shabl=zeros(n(i),2,2);
    for j=1:n(i)
        shabl(j,:,:)=eye(2,2);
    end

    res{i}=shabl;
    res{i}(:,2,1)=ni.*(0:n(i)-1);
    ni=ni*n(i);
end

res{d}=ones(n(d),2);
res{d}(:,2)=ni.*(0:n(d)-1);

res{1}=ones(n(1),2);
res{1}(:,1)=0:n(1)-1;

if (d==1)
    res{1}=res{1}*[1;0];
end;

res = tt_tensor(res); % Bydlocode
return
end
